package analysis;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;
import java.util.regex.Pattern;

/**
 * Analyzes experimental datasets to calculate marker-4 protein abundance ratios.
 * Uses razor theory to select the minimal genome set, then calculates marker-4
 * protein abundances and ratios for each selected genome.
 */
public class DatasetMarkerAnalysis {
	
	// Regex patterns for the 4-category markers (case-insensitive)
	private static final Pattern ANCHOR_PATTERN = Pattern.compile("(?i)\\bribosomal\\s+protein\\b");
	private static final Pattern HIGH_RES_PATTERN = Pattern.compile("(?i)\\b(gyrase|gyr[ab])\\b|\\brna\\s+polymerase\\b|\\bdna\\s+polymerase\\b|\\brpo[abc]\\b|\\bdna[eq]\\b");
	private static final Pattern GROWTH_PATTERN = Pattern.compile("(?i)\\belongation\\s+factor\\b|\\bef[-\\s]?(tu|g|ts|p|4|1a|2)\\b");
	private static final Pattern STABILITY_PATTERN = Pattern.compile("(?i)(aminoacyl[-\\s]?trna|aa[-\\s]?trna|trna)[-\\s]*[^a-zA-Z0-9]{0,2}\\s*(ligase|synthetase)\\b");
	
	// Set of genomes to exclude from results (no function annotation)
	private static final Set<String> EXCLUDED_GENOMES = new HashSet<>(Arrays.asList("PaD", "137", "F0", "F2", "M13", "P22", "ES18", "CRAP"));
	
	// Label to species mapping
	private static final Map<String, String> LABEL_TO_SPECIES = new HashMap<>();
	static {
		LABEL_TO_SPECIES.put("ATN", "Agrobacterium tumefaciens");
		LABEL_TO_SPECIES.put("Am2", "Alteromonas macleodii");
		LABEL_TO_SPECIES.put("BS", "Bacillus subtilis");
		LABEL_TO_SPECIES.put("BXL", "Burkholderia xenovorans");
		LABEL_TO_SPECIES.put("CRH", "Chlamydomonas reinhardtii");
		LABEL_TO_SPECIES.put("CV", "Chromobacterium violaceum");
		LABEL_TO_SPECIES.put("Cup", "Cupriavidus metallidurans");
		LABEL_TO_SPECIES.put("DVH", "Desulfovibrio vulgaris");
		LABEL_TO_SPECIES.put("K12", "Escherichia coli");
		LABEL_TO_SPECIES.put("Ne1", "Nitrosomonas europaeae");
		LABEL_TO_SPECIES.put("Nu1", "Nitrosomonas ureae");
		LABEL_TO_SPECIES.put("NV", "Nitrososphaera viennensis");
		LABEL_TO_SPECIES.put("Nm1", "Nitrosospira multiformis");
		LABEL_TO_SPECIES.put("PaD", "Paracoccus denitrificans");
		LABEL_TO_SPECIES.put("ES18", "Phage ES18");
		LABEL_TO_SPECIES.put("F0", "Phage F0");
		LABEL_TO_SPECIES.put("F2", "Phage F2");
		LABEL_TO_SPECIES.put("M13", "Phage M13");
		LABEL_TO_SPECIES.put("P22", "Phage P22");
		LABEL_TO_SPECIES.put("PD", "Pseudomonas denitrificans");
		LABEL_TO_SPECIES.put("Pfl", "Pseudomonas fluorescens");
		LABEL_TO_SPECIES.put("KF7", "Pseudomonas pseudoalcaligenes");
		LABEL_TO_SPECIES.put("841", "Rhizobium leguminosarum bv. viciae 3841");
		LABEL_TO_SPECIES.put("VF", "Rhizobium leguminosarum bv. viciae VF39");
		LABEL_TO_SPECIES.put("AK199", "Roseobacter sp. AK199");
		LABEL_TO_SPECIES.put("LT2", "Salmonella enterica typhimurium (3 strains combined)");
		LABEL_TO_SPECIES.put("137", "Staphylococcus aureus ATCC 13709");
		LABEL_TO_SPECIES.put("259", "Staphylococcus aureus ATCC 25923");
		LABEL_TO_SPECIES.put("SMS", "Stenotrophomonas maltophilia");
		LABEL_TO_SPECIES.put("HB2", "Thermus Thermophilus");
	}
	
	/**
	 * Extract genome ID from protein ID (e.g., "MGYG000000001_00001" -> "MGYG000000001")
	 */
	private String extractGenomeId(String proteinId) {
		int underscore = proteinId.indexOf('_');
		if (underscore > 0) {
			return proteinId.substring(0, underscore);
		}
		return proteinId;
	}
	
	/**
	 * Load protein-to-category mapping from .faa files
	 */
	private Map<String, String> loadMarkerProteinCategories(File faaFolder) throws IOException {
		Map<String, String> map = new HashMap<>();
		File[] files = faaFolder.listFiles();
		if (files == null) {
			return map;
		}
		for (File f : files) {
			if (!(f.isFile() && (f.getName().endsWith(".faa") || f.getName().endsWith(".fasta")))) {
				continue;
			}
			BufferedReader reader = new BufferedReader(new FileReader(f));
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
				String lineLower = line.toLowerCase();
				String category = categorize(lineLower);
				if (category != null) {
					map.put(proteinId, category);
				}
			}
			reader.close();
		}
		return map;
	}

	private Map<String, String> loadMarkerProteinCategoriesFromFasta(File fastaFile) throws IOException {
		Map<String, String> map = new HashMap<>();
		if (fastaFile == null || !fastaFile.isFile()) {
			return map;
		}
		try (BufferedReader reader = new BufferedReader(new FileReader(fastaFile))) {
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
				String lineLower = line.toLowerCase();
				String category = categorize(lineLower);
				if (category != null) {
					map.put(proteinId, category);
				}
			}
		}
		return map;
	}
	
	private String categorize(String lineLower) {
		if (ANCHOR_PATTERN.matcher(lineLower).find()) return "Anchor";
		if (HIGH_RES_PATTERN.matcher(lineLower).find()) return "HighRes";
		if (GROWTH_PATTERN.matcher(lineLower).find()) return "Growth";
		if (STABILITY_PATTERN.matcher(lineLower).find()) return "Housekeeping";
		return null;
	}
	
	/**
	 * Class to hold extracted protein information
	 */
	public static class ProteinInfo {
		public String reference;
		public String taxa;
		public String markerCategory; // "Anchor", "HighRes", "Growth", "Housekeeping", or "None"
		
		public ProteinInfo(String reference, String taxa, String markerCategory) {
			this.reference = reference;
			this.taxa = taxa;
			this.markerCategory = markerCategory;
		}
		
		@Override
		public String toString() {
			return String.format("Reference: %s, Taxa: %s, Marker: %s", reference, taxa, markerCategory);
		}
	}
	
	/**
	 * Extract protein information from a FASTA header line
	 * @param headerLine The FASTA header line (starting with ">")
	 * @return ProteinInfo object containing reference, taxa, and marker category
	 */
	public ProteinInfo extractProteinInfo(String headerLine) {
		if (headerLine == null || !headerLine.startsWith(">")) {
			return null;
		}
		
		String trimmed = headerLine.substring(1).trim();
		if (trimmed.isEmpty()) {
			return null;
		}
		
		// Extract reference (first token)
		String[] parts = trimmed.split("\\s+", 2);
		String reference = parts[0];
		
		// Extract label (first part before underscore)
		String taxa = null;
		int underscoreIndex = reference.indexOf('_');
		if (underscoreIndex > 0) {
			String label = reference.substring(0, underscoreIndex);
			taxa = LABEL_TO_SPECIES.get(label);
		} else {
			// No underscore, check if entire reference is a label
			taxa = LABEL_TO_SPECIES.get(reference);
		}
		
		// If label not found in mapping, try to extract from header
		if (taxa == null || taxa.isEmpty()) {
			taxa = extractTaxa(trimmed);
		}
		
		// Normalize taxa name (e.g., combine Stenotrophomonas variants)
		if (taxa != null && !taxa.isEmpty()) {
			taxa = normalizeTaxa(taxa);
		}
		
		// Extract marker category
		String lineLower = headerLine.toLowerCase();
		String category = categorize(lineLower);
		String markerCategory = (category != null) ? category : "None";
		
		return new ProteinInfo(reference, taxa, markerCategory);
	}
	
	/**
	 * Normalize taxa names to combine variants
	 * @param taxa Original taxa name
	 * @return Normalized taxa name
	 */
	private String normalizeTaxa(String taxa) {
		if (taxa == null || taxa.isEmpty()) {
			return taxa;
		}
		
		// Combine "Stenotrophomonas maltophilia" to "Stenotrophomonas"
		// (SMS label maps to "Stenotrophomonas maltophilia" which will be normalized here)
		if (taxa.equals("Stenotrophomonas maltophilia")) {
			return "Stenotrophomonas";
		}
		
		return taxa;
	}
	
	/**
	 * Extract taxa information from FASTA header
	 * Supports two formats:
	 * 1. OS=Taxa name (strain info) - extract text between OS= and (
	 * 2. [Taxa name] - extract text between brackets
	 */
	private String extractTaxa(String header) {
		// Try format 1: OS=Taxa name (strain info)
		int osIndex = header.indexOf("OS=");
		if (osIndex >= 0) {
			int start = osIndex + 3; // Skip "OS="
			int end = header.indexOf("(", start);
			if (end < 0) {
				// No opening parenthesis, try to find end of line or next space-separated field
				end = header.length();
				// Look for common field separators
				int nextField = header.indexOf(" GN=", start);
				if (nextField > 0 && nextField < end) end = nextField;
				nextField = header.indexOf(" PE=", start);
				if (nextField > 0 && nextField < end) end = nextField;
				nextField = header.indexOf(" SV=", start);
				if (nextField > 0 && nextField < end) end = nextField;
			}
			if (end > start) {
				String taxa = header.substring(start, end).trim();
				if (!taxa.isEmpty()) {
					return taxa;
				}
			}
		}
		
		// Try format 2: [Taxa name]
		int bracketStart = header.lastIndexOf('[');
		int bracketEnd = header.indexOf(']', bracketStart);
		if (bracketStart >= 0 && bracketEnd > bracketStart) {
			String taxa = header.substring(bracketStart + 1, bracketEnd).trim();
			if (!taxa.isEmpty()) {
				return taxa;
			}
		}
		
		// If no taxa found, return empty string
		return "";
	}
	
	/**
	 * Extract marker proteins from a FASTA database to create a reduced database
	 * Only includes marker proteins (Anchor, HighRes, Growth, Housekeeping) and excludes proteins from EXCLUDED_GENOMES
	 * @param fastaFile Input FASTA file
	 * @param outputFile Output FASTA file containing only marker proteins
	 */
	public void extractMarkerProteinsFromFasta(File fastaFile, File outputFile) throws IOException {
		if (!fastaFile.exists() || !fastaFile.isFile()) {
			throw new IOException("FASTA file does not exist: " + fastaFile.getAbsolutePath());
		}
		
		System.out.println("=== Extracting Marker Proteins from FASTA Database ===");
		System.out.println("Input file: " + fastaFile.getPath());
		System.out.println("Output file: " + outputFile.getPath());
		System.out.println("Excluding genomes: " + EXCLUDED_GENOMES);
		
		PrintWriter writer = new PrintWriter(outputFile);
		int totalProteins = 0;
		int markerProteins = 0;
		int excludedProteins = 0;
		Map<String, Integer> categoryCounts = new HashMap<>();
		
		try (BufferedReader reader = new BufferedReader(new FileReader(fastaFile))) {
			String line;
			String currentHeader = null;
			StringBuilder currentSequence = new StringBuilder();
			boolean inMarkerProtein = false;
			boolean excludeThisProtein = false;
			
			while ((line = reader.readLine()) != null) {
				if (line.startsWith(">")) {
					// Process previous protein if any
					if (currentHeader != null && inMarkerProtein && !excludeThisProtein && currentSequence.length() > 0) {
						// Write the marker protein to output
						writer.println(currentHeader);
						// Write sequence in 60-character lines (FASTA format)
						String sequence = currentSequence.toString();
						for (int i = 0; i < sequence.length(); i += 60) {
							int end = Math.min(i + 60, sequence.length());
							writer.println(sequence.substring(i, end));
						}
						markerProteins++;
					}
					
					// Start new protein
					currentHeader = line;
					currentSequence = new StringBuilder();
					
					// Skip if title contains CRAP
					String lineUpper = line.toUpperCase();
					if (lineUpper.contains("CRAP_")) {
						excludedProteins++;
						totalProteins++;
						excludeThisProtein = true;
						inMarkerProtein = false;
						continue;
					}
					
					// Extract protein info
					ProteinInfo info = extractProteinInfo(line);
					if (info != null) {
						totalProteins++;
						
						// Check if it's a marker protein (category != "None")
						inMarkerProtein = !info.markerCategory.equals("None");
						
						// Check if genome should be excluded
						// Extract identifier (first part before underscore) as genome identifier
						String genome = null;
						int underscoreIndex = info.reference.indexOf('_');
						if (underscoreIndex > 0) {
							genome = info.reference.substring(0, underscoreIndex);
						} else {
							genome = info.reference;
						}
						
						excludeThisProtein = EXCLUDED_GENOMES.contains(genome);
						
						if (excludeThisProtein) {
							excludedProteins++;
						} else if (inMarkerProtein) {
							categoryCounts.merge(info.markerCategory, 1, Integer::sum);
						}
					} else {
						excludeThisProtein = true;
						inMarkerProtein = false;
					}
				} else {
					// Sequence line
					if (currentHeader != null) {
						currentSequence.append(line.trim());
					}
				}
			}
			
			// Process last protein
			if (currentHeader != null && inMarkerProtein && !excludeThisProtein && currentSequence.length() > 0) {
				writer.println(currentHeader);
				String sequence = currentSequence.toString();
				for (int i = 0; i < sequence.length(); i += 60) {
					int end = Math.min(i + 60, sequence.length());
					writer.println(sequence.substring(i, end));
				}
				markerProteins++;
			}
		}
		
		writer.close();
		
		System.out.println("\n=== Extraction Summary ===");
		System.out.println("Total proteins processed: " + totalProteins);
		System.out.println("Marker proteins extracted: " + markerProteins);
		System.out.println("Proteins excluded (from EXCLUDED_GENOMES or CRAP): " + excludedProteins);
		System.out.println("\nMarker protein counts by category:");
		for (Map.Entry<String, Integer> entry : categoryCounts.entrySet()) {
			System.out.println("  " + entry.getKey() + ": " + entry.getValue());
		}
		System.out.println("\nReduced database written to: " + outputFile.getPath());
	}
	
	/**
	 * Extract protein information from a FASTA file and write to output file
	 * @param fastaFile Input FASTA file
	 * @param outputFile Output TSV file with columns: Reference, Taxa, Marker_Category
	 */
	public void extractProteinInfoFromFasta(File fastaFile, File outputFile) throws IOException {
		if (!fastaFile.exists() || !fastaFile.isFile()) {
			throw new IOException("FASTA file does not exist: " + fastaFile.getAbsolutePath());
		}
		
		// Statistics tracking: taxa -> marker category -> count
		Map<String, Map<String, Integer>> taxaStats = new HashMap<>();
		Set<String> allTaxa = new HashSet<>();
		int totalProteins = 0;
		
		// Collect examples of "Unknown" taxa headers (max 20 examples)
		ArrayList<String> unknownTaxaExamples = new ArrayList<>();
		int maxExamples = 200;
		
		PrintWriter writer = new PrintWriter(outputFile);
		writer.println("Reference\tTaxa\tMarker_Category");
		
		HashSet<String> startSet = new HashSet<String>();
		try (BufferedReader reader = new BufferedReader(new FileReader(fastaFile))) {
			String line;
			int skippedCount = 0;
			while ((line = reader.readLine()) != null) {
				if (line.startsWith(">")) {
					// Skip if title contains BOVIN or HUMAN
					String lineUpper = line.toUpperCase();
					if (lineUpper.contains("CRAP_")) {
						skippedCount++;
						continue;
					}
					
					String starter = line.substring(1, line.indexOf("_"));
					startSet.add(starter);
					
					ProteinInfo info = extractProteinInfo(line);
					if (info != null) {
						String taxa = info.taxa.isEmpty() ? "Unknown" : info.taxa;
						String marker = info.markerCategory;
						
						writer.printf("%s\t%s\t%s\n", info.reference, taxa, marker);
						
						// Collect examples of Unknown taxa
						if (taxa.equals("Unknown") && unknownTaxaExamples.size() < maxExamples) {
							unknownTaxaExamples.add(line);
						}
						
						// Update statistics
						allTaxa.add(taxa);
						totalProteins++;
						taxaStats.computeIfAbsent(taxa, k -> new HashMap<>())
							.merge(marker, 1, Integer::sum);
					}
				}
			}
			
			if (skippedCount > 0) {
				System.out.println("Skipped " + skippedCount + " proteins containing BOVIN or HUMAN");
			}
		}
		
		System.out.println(startSet.size());
		for(String ss: startSet) {
			System.out.println(ss);
		}
		
		writer.close();
		System.out.println("Extracted protein information from " + fastaFile.getName() + 
			" to " + outputFile.getName());
		System.out.println("Total proteins processed: " + totalProteins);
		System.out.println("Total unique taxa: " + allTaxa.size());
		
		// Print examples of Unknown taxa
		int unknownCount = taxaStats.getOrDefault("Unknown", new HashMap<>()).values().stream()
			.mapToInt(Integer::intValue).sum();
		if (unknownCount > 0) {
			System.out.println("\n=== Unknown Taxa Examples ===");
			System.out.println("Total proteins with Unknown taxa: " + unknownCount);
			System.out.println("Showing up to " + maxExamples + " examples:");
			System.out.println();
			for (int i = 0; i < unknownTaxaExamples.size(); i++) {
				System.out.println("Example " + (i + 1) + ":");
				System.out.println("  " + unknownTaxaExamples.get(i));
				System.out.println();
			}
			if (unknownCount > maxExamples) {
				System.out.println("... and " + (unknownCount - maxExamples) + " more Unknown taxa entries");
			}
		}
		
		// Write statistics file
		File statsFile = new File(outputFile.getParent(), 
			outputFile.getName().replaceFirst("\\.[^.]+$", "_statistics.txt"));
		writeTaxaStatistics(taxaStats, allTaxa, statsFile);
		System.out.println("Statistics written to: " + statsFile.getName());
	}
	
	/**
	 * Write taxa statistics: total taxa count and marker protein fractions per taxa
	 * @param taxaStats Map of taxa -> marker category -> count
	 * @param allTaxa Set of all taxa
	 * @param statsFile Output file for statistics
	 */
	private void writeTaxaStatistics(Map<String, Map<String, Integer>> taxaStats, 
			Set<String> allTaxa, File statsFile) throws IOException {
		
		PrintWriter writer = new PrintWriter(statsFile);
		writer.println("=== Protein Database Statistics ===");
		writer.println();
		writer.println("Total number of taxa: " + allTaxa.size());
		writer.println();
		writer.println("=== Marker Protein Fractions by Taxa ===");
		writer.println();
		
		// Header
		writer.println("Taxa\tTotal_Proteins\tAnchor\tHighRes\tGrowth\tHousekeeping\tNone\t" +
			"Anchor_Fraction\tHighRes_Fraction\tGrowth_Fraction\tHousekeeping_Fraction\tNone_Fraction");
		
		// Sort taxa for consistent output
		ArrayList<String> sortedTaxa = new ArrayList<>(allTaxa);
		Collections.sort(sortedTaxa);
		
		for (String taxa : sortedTaxa) {
			Map<String, Integer> markerCounts = taxaStats.get(taxa);
			if (markerCounts == null) {
				continue;
			}
			
			// Get counts for each category
			int anchorCount = markerCounts.getOrDefault("Anchor", 0);
			int highResCount = markerCounts.getOrDefault("HighRes", 0);
			int growthCount = markerCounts.getOrDefault("Growth", 0);
			int housekeepingCount = markerCounts.getOrDefault("Housekeeping", 0);
			int noneCount = markerCounts.getOrDefault("None", 0);
			int totalCount = anchorCount + highResCount + growthCount + housekeepingCount + noneCount;
			
			// Calculate fractions
			double anchorFraction = totalCount > 0 ? (double) anchorCount / totalCount : 0.0;
			double highResFraction = totalCount > 0 ? (double) highResCount / totalCount : 0.0;
			double growthFraction = totalCount > 0 ? (double) growthCount / totalCount : 0.0;
			double housekeepingFraction = totalCount > 0 ? (double) housekeepingCount / totalCount : 0.0;
			double noneFraction = totalCount > 0 ? (double) noneCount / totalCount : 0.0;
			
			writer.printf("%s\t%d\t%d\t%d\t%d\t%d\t%d\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\n",
				taxa, totalCount, anchorCount, highResCount, growthCount, housekeepingCount, noneCount,
				anchorFraction, highResFraction, growthFraction, housekeepingFraction, noneFraction);
		}
		
		writer.close();
	}
	
	/**
	 * Statistics tracking class
	 */
	private static class AnalysisStatistics {
		int totalPeptidesInFile = 0;
		int peptidesWithMGYG = 0;
		int peptidesWithIntensities = 0;
		int peptidesMappedToGenomes = 0;
		int peptidesUsedInRazor = 0;
		int peptidesUsedInAbundance = 0;
		int distinctPeptidesUsedInAbundance = 0;
		int totalGenomesFound = 0;
		int genomesSelectedByRazor = 0;
		int decoyProteinsSkipped = 0;
		int hostProteinsSkipped = 0;
		int otherProteinsSkipped = 0;
		
		// Detectability statistics
		int uniquePeptides = 0;  // Peptides in 1 genome
		int lowSpecificityPeptides = 0;  // Peptides in 2-5 genomes
		int mediumSpecificityPeptides = 0;  // Peptides in 6-20 genomes
		int highSpecificityPeptides = 0;  // Peptides in >20 genomes
		
		void writeToFile(File outputFile, Map<String, Set<String>> peptideToGenomes) throws IOException {
			// Calculate detectability statistics
			uniquePeptides = 0;
			lowSpecificityPeptides = 0;
			mediumSpecificityPeptides = 0;
			highSpecificityPeptides = 0;
			
			for (Set<String> genomes : peptideToGenomes.values()) {
				int count = genomes.size();
				if (count == 1) {
					uniquePeptides++;
				} else if (count <= 5) {
					lowSpecificityPeptides++;
				} else if (count <= 20) {
					mediumSpecificityPeptides++;
				} else {
					highSpecificityPeptides++;
				}
			}
			
			PrintWriter writer = new PrintWriter(outputFile);
			writer.println("=== Peptide Mapping Statistics ===");
			writer.println();
			writer.println("Input Statistics:");
			writer.println("  Total peptides in file: " + totalPeptidesInFile);
			writer.println("  Peptides with MGYG proteins: " + peptidesWithMGYG);
			writer.println("  Peptides with valid intensities: " + peptidesWithIntensities);
			writer.println();
			writer.println("Protein Filtering:");
			writer.println("  Decoy proteins skipped (REV): " + decoyProteinsSkipped);
			writer.println("  Host proteins skipped (sp): " + hostProteinsSkipped);
			writer.println("  Other proteins skipped: " + otherProteinsSkipped);
			writer.println();
			writer.println("Genome Mapping:");
			writer.println("  Peptides mapped to genomes: " + peptidesMappedToGenomes);
			writer.println("  Total genomes found: " + totalGenomesFound);
			writer.println();
			writer.println("Peptide Detectability:");
			writer.println("  Unique peptides (1 genome): " + uniquePeptides);
			writer.println("  Low specificity (2-5 genomes): " + lowSpecificityPeptides);
			writer.println("  Medium specificity (6-20 genomes): " + mediumSpecificityPeptides);
			writer.println("  High specificity (>20 genomes): " + highSpecificityPeptides);
			double uniquePct = peptidesMappedToGenomes > 0 ? 
				(100.0 * uniquePeptides / peptidesMappedToGenomes) : 0.0;
			writer.printf("  Unique peptide percentage: %.2f%%\n", uniquePct);
			writer.println();
			writer.println("Razor Selection:");
			writer.println("  Peptides used in razor: " + peptidesUsedInRazor);
			writer.println("  Genomes selected by razor: " + genomesSelectedByRazor);
			writer.println();
			writer.println("Abundance Calculation:");
			writer.println("  Peptides used in abundance calculation (all): " + peptidesUsedInAbundance);
			writer.println("  Distinct peptides used in abundance calculation: " + distinctPeptidesUsedInAbundance);
			double distinctPct = peptidesUsedInAbundance > 0 ? 
				(100.0 * distinctPeptidesUsedInAbundance / peptidesUsedInAbundance) : 0.0;
			writer.printf("  Distinct peptide percentage: %.2f%%\n", distinctPct);
			writer.println();
			writer.println("Coverage:");
			double coverage = totalPeptidesInFile > 0 ? 
				(100.0 * peptidesUsedInAbundance / totalPeptidesInFile) : 0.0;
			writer.printf("  Peptide coverage: %.2f%% (%d / %d)\n", 
				coverage, peptidesUsedInAbundance, totalPeptidesInFile);
			writer.close();
		}
	}
	
	/**
	 * Strip modifications from peptide sequence (e.g., "(Acetyl)K..." -> "K...")
	 * Handles nested parentheses like "(Succinyl_13C(4))" by removing balanced parentheses
	 */
	private String stripModifications(String sequence) {
		if (sequence == null || sequence.isEmpty()) {
			return sequence;
		}
		
		StringBuilder result = new StringBuilder();
		int i = 0;
		int len = sequence.length();
		
		while (i < len) {
			char c = sequence.charAt(i);
			
			if (c == '(') {
				// Found opening parenthesis - find matching closing parenthesis
				int depth = 1;
				int j = i + 1;
				
				// Find the matching closing parenthesis by counting depth
				while (j < len && depth > 0) {
					char next = sequence.charAt(j);
					if (next == '(') {
						depth++;
					} else if (next == ')') {
						depth--;
					}
					j++;
				}
				
				// Skip the entire modification (from '(' to matching ')')
				i = j;
			} else {
				// Regular character - keep it
				result.append(c);
				i++;
			}
		}
		
		return result.toString().trim();
	}
	
	/**
	 * Read peptide-intensity file and extract peptides with their intensities and protein mappings
	 * Returns: Map<peptide, PeptideData>
	 */
	private static class PeptideData {
		String originalSequence;
		String strippedSequence;
		Set<String> proteinIds;
		Map<String, Double> intensities;
		
		PeptideData(String original, String stripped, Set<String> proteins, Map<String, Double> intensities) {
			this.originalSequence = original;
			this.strippedSequence = stripped;
			this.proteinIds = proteins;
			this.intensities = intensities;
		}
	}
	
	/**
	 * Read peptide-intensity file and extract peptides with their intensities and protein mappings
	 * @param filterMGYGOnly If true, only keep proteins starting with MGYG. If false, keep all proteins except REV decoys and sp host.
	 */
	private Map<String, PeptideData> readPeptideIntensities(File peptideFile, AnalysisStatistics stats, boolean filterMGYGOnly) throws IOException {
		Map<String, PeptideData> peptideData = new HashMap<>();
		BufferedReader reader = new BufferedReader(new FileReader(peptideFile));
		
		// Read header
		String headerLine = reader.readLine();
		if (headerLine == null) {
			reader.close();
			return peptideData;
		}
		
		String[] headers = headerLine.split("\t");
		ArrayList<String> intensityColumns = new ArrayList<>();
		int sequenceIdx = -1;
		int proteinsIdx = -1;
		
		// Find relevant column indices
		for (int i = 0; i < headers.length; i++) {
			String h = headers[i].trim();
			if (h.equalsIgnoreCase("Sequence") || h.equalsIgnoreCase("peptide")) {
				sequenceIdx = i;
			} else if (h.equalsIgnoreCase("Proteins")) {
				proteinsIdx = i;
			} else if (h.toLowerCase().contains("intensity") && !h.equalsIgnoreCase("Intensity")) {
				// Skip the total "Intensity" column, get sample-specific ones
				intensityColumns.add(h);
			} else if (h.toLowerCase().contains(".mzml")) {
				// Skip the total "Intensity" column, get sample-specific ones
				intensityColumns.add(h);
			}
		}
		
		if (sequenceIdx < 0 || proteinsIdx < 0) {
			System.err.println("Warning: Could not find Sequence or Proteins column");
			reader.close();
			return peptideData;
		}
		
		// Read data
		String line;
		while ((line = reader.readLine()) != null) {
			stats.totalPeptidesInFile++;
			
			String[] parts = line.split("\t");
			if (parts.length <= Math.max(sequenceIdx, proteinsIdx)) {
				continue;
			}
			
			String sequence = parts[sequenceIdx].trim();
			String proteins = parts[proteinsIdx].trim();
			
			if (sequence.isEmpty() || proteins.isEmpty()) {
				continue;
			}
			
			// Strip modifications for matching
			String strippedSequence = stripModifications(sequence);
			
			// Parse protein IDs (semicolon-separated) and filter
			Set<String> proteinIds = new HashSet<>();
			Set<String> allProteinIds = new HashSet<>(); // Track all proteins for statistics
			for (String proteinId : proteins.split(";")) {
				proteinId = proteinId.trim();
				if (!proteinId.isEmpty()) {
					allProteinIds.add(proteinId);
					
					// Always exclude REV decoys and sp host proteins
					if (proteinId.startsWith("REV")) {
						stats.decoyProteinsSkipped++;
						continue;
					} else if (proteinId.startsWith("sp")) {
						stats.hostProteinsSkipped++;
						continue;
					}
					
					// Filter by MGYG if requested, otherwise accept all non-decoys
					if (filterMGYGOnly) {
						if (proteinId.startsWith("MGYG")) {
							proteinIds.add(proteinId);
						} else {
							stats.otherProteinsSkipped++;
						}
					} else {
						// Accept all proteins (already filtered out REV and sp)
						proteinIds.add(proteinId);
					}
				}
			}
			
			if (!proteinIds.isEmpty()) {
				if (filterMGYGOnly) {
					stats.peptidesWithMGYG++;
				} else {
					stats.peptidesWithMGYG++; // Reuse this counter for "peptides with valid proteins"
				}
			}
			
			// Extract intensities for each sample
			Map<String, Double> intensities = new HashMap<>();
			for (String intensityCol : intensityColumns) {
				// Find column index for this intensity column
				int colIdx = -1;
				for (int i = 0; i < headers.length; i++) {
					if (headers[i].equals(intensityCol)) {
						colIdx = i;
						break;
					}
				}
				if (colIdx >= 0 && colIdx < parts.length) {
					try {
						double intensity = Double.parseDouble(parts[colIdx].trim());
						if (intensity > 0) {
							intensities.put(intensityCol, intensity);
						}
					} catch (NumberFormatException e) {
						// Skip invalid numbers
					}
				}
			}
			
			// Only store peptides that have valid proteins and intensities
			if (!intensities.isEmpty() && !proteinIds.isEmpty()) {
				stats.peptidesWithIntensities++;
				// Store using stripped sequence as key (for matching with predicted files)
				// If strippedSequence already exists, merge the data (combine intensities and protein IDs)
				if (peptideData.containsKey(strippedSequence)) {
					// Merge with existing entry
					PeptideData existing = peptideData.get(strippedSequence);
					// Merge protein IDs
					existing.proteinIds.addAll(proteinIds);
					// Combine intensities (add values for each sample)
					for (Map.Entry<String, Double> intensityEntry : intensities.entrySet()) {
						String sample = intensityEntry.getKey();
						double newIntensity = intensityEntry.getValue();
						existing.intensities.merge(sample, newIntensity, Double::sum);
					}
				} else {
					// Create new entry
					peptideData.put(strippedSequence, new PeptideData(sequence, strippedSequence, proteinIds, intensities));
				}
			}
		}
		
		reader.close();
		return peptideData;
	}
	
	/**
	 * Map peptides to genomes from peptide data (using protein IDs)
	 * Returns: Map<peptide, Set<genomeId>>
	 */
	private Map<String, Set<String>> mapPeptidesToGenomes(Map<String, PeptideData> peptideData) {
		Map<String, Set<String>> peptideToGenomes = new HashMap<>();
		
		for (Map.Entry<String, PeptideData> entry : peptideData.entrySet()) {
			String peptide = entry.getKey();
			PeptideData data = entry.getValue();
			
			Set<String> genomes = new HashSet<>();
			for (String proteinId : data.proteinIds) {
				String genomeId = extractGenomeId(proteinId);
				genomes.add(genomeId);
			}
			
			if (!genomes.isEmpty()) {
				peptideToGenomes.put(peptide, genomes);
			}
		}
		
		return peptideToGenomes;
	}
	
	/**
	 * Mapping result class for predicted file mapping
	 */
	private static class PredictedMappingResult {
		Map<String, Set<String>> peptideToGenomes;
		Map<String, String> proteinToGenome;  // Maps protein ID to genome ID
		Map<String, Double> peptideToDetectability;  // Maps peptide to detectability (from predicted files)
		
		PredictedMappingResult() {
			this.peptideToGenomes = new HashMap<>();
			this.proteinToGenome = new HashMap<>();
			this.peptideToDetectability = new HashMap<>();
		}
	}
	
	/**
	 * Map peptides to genomes by reading from predicted files folder
	 * Extracts genome ID from filename (e.g., "MGYG000000005.predicted.tsv" -> "MGYG000000005")
	 * Only maps peptides that are in the peptideData set (for efficiency)
	 * Also builds protein-to-genome mapping for abundance calculations
	 * Returns: PredictedMappingResult with both peptide-to-genomes and protein-to-genome mappings
	 */
	private PredictedMappingResult mapPeptidesToGenomesFromPredictedFiles(
			File predictedFolder, Map<String, PeptideData> peptideData) throws IOException {
		
		PredictedMappingResult result = new PredictedMappingResult();
		
		// First, create a set of all stripped peptide sequences we're looking for (for fast lookup)
		Set<String> targetPeptides = new HashSet<>(peptideData.keySet());
		
		// Get all .predicted.tsv files in the folder
		File[] files = predictedFolder.listFiles();
		if (files == null) {
			System.err.println("Warning: Predicted folder is empty or cannot be read: " + predictedFolder.getPath());
			return result;
		}
		
		System.out.println("Reading predicted files from: " + predictedFolder.getPath());
		int processedFiles = 0;
		
		for (File file : files) {
			if (!file.isFile() || !file.getName().endsWith(".predicted.tsv")) {
				continue;
			}
			
			// Extract genome ID from filename (e.g., "MGYG000000005.predicted.tsv" -> "MGYG000000005")
			String fileName = file.getName();
			String genomeId = fileName.substring(0, fileName.length() - ".predicted.tsv".length());
			
			try (BufferedReader reader = new BufferedReader(new FileReader(file))) {
				// Read header
				String headerLine = reader.readLine();
				if (headerLine == null) {
					continue;
				}
				
				String[] headers = headerLine.split("\t");
				int peptideSeqIdx = -1;
				int proteinIdIdx = -1;
				int detectabilityIdx = -1;
				
				// Find column indices
				for (int i = 0; i < headers.length; i++) {
					String h = headers[i].trim();
					if (h.equals("Peptide sequence")) {
						peptideSeqIdx = i;
					} else if (h.equals("Protein id")) {
						proteinIdIdx = i;
					} else if (h.equals("Peptide detectability")) {
						detectabilityIdx = i;
					}
				}
				
				if (peptideSeqIdx < 0 || proteinIdIdx < 0) {
					System.err.println("Warning: Could not find required columns in " + fileName);
					continue;
				}
				
				// Read data lines
				String line;
				while ((line = reader.readLine()) != null) {
					String[] parts = line.split("\t");
					if (parts.length <= Math.max(peptideSeqIdx, proteinIdIdx)) {
						continue;
					}
					
					String peptideSequence = parts[peptideSeqIdx].trim();
					String proteinId = parts[proteinIdIdx].trim();
					
					// Filter: only include proteins starting with MGYG (exclude REV decoys and sp host)
					if (!proteinId.startsWith("MGYG")) {
						continue;
					}
					
					// Build protein-to-genome mapping (for abundance calculations)
					result.proteinToGenome.put(proteinId, genomeId);
					
					// Only process peptides that are in our target set
					if (targetPeptides.contains(peptideSequence)) {
						result.peptideToGenomes.computeIfAbsent(peptideSequence, k -> new HashSet<>()).add(genomeId);
						
						// Read detectability value if column exists
						if (detectabilityIdx >= 0 && detectabilityIdx < parts.length) {
							try {
								double detectability = Double.parseDouble(parts[detectabilityIdx].trim());
								// Keep the maximum detectability for each peptide (in case it appears multiple times)
								Double existingDetectability = result.peptideToDetectability.get(peptideSequence);
								if (existingDetectability == null || detectability > existingDetectability) {
									result.peptideToDetectability.put(peptideSequence, detectability);
								}
							} catch (NumberFormatException e) {
								// Skip invalid detectability values
							}
						}
					}
				}
			} catch (IOException e) {
				System.err.println("Error reading file " + fileName + ": " + e.getMessage());
			}
			
			processedFiles++;
			if (processedFiles % 100 == 0) {
				System.out.println("  Processed " + processedFiles + " predicted files...");
			}
		}
		
		System.out.println("Finished processing " + processedFiles + " predicted files");
		System.out.println("Mapped " + result.peptideToGenomes.size() + " peptides to genomes");
		System.out.println("Mapped " + result.proteinToGenome.size() + " proteins to genomes");
		System.out.println("Found detectability values for " + result.peptideToDetectability.size() + " peptides");
		
		return result;
	}
	
	/**
	 * Apply Occam's Razor to select minimal genome set
	 * Based on MetaPilotTask.applyOccamsRazor
	 */
	private HashSet<String> applyOccamsRazor(HashMap<String, HashSet<String>> genomeToPeptides,
			Set<String> totalPeptides, double threshold, double stepThreshold, int maximum, File out) {
		
		int totalPepCount = totalPeptides.size();
		String[] genomes = genomeToPeptides.keySet().toArray(new String[genomeToPeptides.size()]);
		
		// Initial sort by peptide count (descending)
		Arrays.sort(genomes, (g1, g2) -> {
			int s1 = genomeToPeptides.get(g1).size();
			int s2 = genomeToPeptides.get(g2).size();
			return s2 - s1;
		});
		
		// Iterative re-sorting based on actual contributions (5 iterations)
		for (int iteratorCount = 0; iteratorCount < 5; iteratorCount++) {
			HashSet<String> totalSet = new HashSet<String>();
			HashMap<String, Integer> orderMap = new HashMap<String, Integer>();
			for (int i = 0; i < genomes.length; i++) {
				int currentSize = totalSet.size();
				totalSet.addAll(genomeToPeptides.get(genomes[i]));
				int increaseSize = totalSet.size() - currentSize;
				orderMap.put(genomes[i], increaseSize);
			}
			
			// Re-sort by actual increment contributions
			Arrays.sort(genomes, (g1, g2) -> {
				return orderMap.get(g2) - orderMap.get(g1);
			});
		}
		
		double[] percentage = new double[genomes.length];
		HashSet<String> selectedGenomes = new HashSet<String>();
		HashSet<String> totalSet = new HashSet<String>();
		int breakPoint1 = -1;
		int breakPoint2 = -1;
		
		// Write ranking statistics if output file is provided
		if (out != null) {
			try (PrintWriter writer = new PrintWriter(out)) {
				writer.println("Rank\tGenome\tPeptide_count\tTotal_covered_peptide_count\tPercentage\tIncrement\tIncrement_ratio");
				for (int i = 0; i < genomes.length; i++) {
					int currentSize = totalSet.size();
					totalSet.addAll(genomeToPeptides.get(genomes[i]));
					int add = totalSet.size() - currentSize;
					double additional = ((double) add / (double) totalPepCount);
					
					percentage[i] = ((double) totalSet.size() / (double) totalPepCount);
					
					writer.println(i + "\t" + genomes[i] + "\t" + genomeToPeptides.get(genomes[i]).size() + "\t"
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
				System.err.println("Error writing ranking file: " + e.getMessage());
			}
		} else {
			// Calculate breakpoints without writing to file
			for (int i = 0; i < genomes.length; i++) {
				int currentSize = totalSet.size();
				totalSet.addAll(genomeToPeptides.get(genomes[i]));
				int add = totalSet.size() - currentSize;
				double additional = ((double) add / (double) totalPepCount);
				
				percentage[i] = ((double) totalSet.size() / (double) totalPepCount);
				
				if (add > 1) {
					breakPoint1 = i;
				}
				if (additional > stepThreshold) {
					breakPoint2 = i;
				}
			}
		}
		
		// Reset totalSet for actual selection
		totalSet.clear();
		
		// Select genomes
		for (int i = 0; i < genomes.length; i++) {
			selectedGenomes.add(genomes[i]);
			
			if (percentage[i] >= threshold) {
				break;
			}
			if (threshold < 1.0) {
				if (i == breakPoint1 || i == breakPoint2) {
					break;
				}
				if (selectedGenomes.size() > maximum) {
					break;
				}
			}
		}
		
		return selectedGenomes;
	}
	
	/**
	 * Calculate marker-4 protein abundances for selected genomes
	 * @param proteinToGenome If provided, uses this mapping instead of extracting genome ID from protein name
	 */
	private void calculateMarkerAbundances(Map<String, PeptideData> peptideData, File faaFolder,
			HashSet<String> selectedGenomes, File outputFile, AnalysisStatistics stats,
			Map<String, String> proteinToGenome) throws IOException {
		
		System.out.println("Loading protein-to-category mapping...");
		Map<String, String> proteinToCategory = loadMarkerProteinCategories(faaFolder);
		System.out.println("Selected genome count: "+selectedGenomes.size());
		// Map: sample -> genome -> category -> abundance
		Map<String, Map<String, Map<String, Double>>> sampleGenomeCategoryAbundance = new HashMap<>();
		// Map: sample -> genome -> total abundance (all proteins)
		Map<String, Map<String, Double>> sampleGenomeTotalAbundance = new HashMap<>();
		
		System.out.println("Calculating abundances...");
		// Calculate abundances
		for (Map.Entry<String, PeptideData> entry : peptideData.entrySet()) {
			String peptide = entry.getKey();
			PeptideData data = entry.getValue();
			
			// Group proteins by genome
			Map<String, Set<String>> genomeProteins = new HashMap<>();
			for (String proteinId : data.proteinIds) {
				String genomeId;
				if (proteinToGenome != null) {
					// Use provided protein-to-genome mapping
					genomeId = proteinToGenome.get(proteinId);
					if (genomeId == null) {
						continue;  // Skip if protein not found in mapping
					}
				} else {
					// Extract genome ID from protein name
					genomeId = extractGenomeId(proteinId);
				}
				
				if (selectedGenomes.contains(genomeId)) {
					genomeProteins.computeIfAbsent(genomeId, k -> new HashSet<>()).add(proteinId);
				}
			}
			
			if (genomeProteins.isEmpty()) {
				continue;
			}
			
			// Track peptides used in abundance calculation
			stats.peptidesUsedInAbundance++;
			
			// For each sample
			for (Map.Entry<String, Double> intensityEntry : data.intensities.entrySet()) {
				String sample = intensityEntry.getKey();
				double intensity = intensityEntry.getValue();
				
				// For each genome that has this peptide
				for (Map.Entry<String, Set<String>> genomeEntry : genomeProteins.entrySet()) {
					String genomeId = genomeEntry.getKey();
					Set<String> proteins = genomeEntry.getValue();
					
					// Distribute intensity equally among proteins in this genome
					// (razor principle: if peptide maps to multiple proteins in same genome, split equally)
					double intensityPerProtein = intensity / proteins.size();
					
					// Add to total genome abundance (use full intensity, not per-protein)
					sampleGenomeTotalAbundance.computeIfAbsent(sample, k -> new HashMap<>())
						.merge(genomeId, intensity, Double::sum);
					
					// Add to category abundances (distribute per protein)
					for (String proteinId : proteins) {
						String category = proteinToCategory.get(proteinId);
						if (category != null) {
							sampleGenomeCategoryAbundance.computeIfAbsent(sample, k -> new HashMap<>())
								.computeIfAbsent(genomeId, k -> new HashMap<>())
								.merge(category, intensityPerProtein, Double::sum);
						}
					}
				}
			}
		}
		
		// Write results
		System.out.println("Writing results...");
		PrintWriter writer = new PrintWriter(outputFile);
		writer.println("Sample\tGenome\tCategory\tAbundance\tTotal_Genome_Abundance\tRatio");
		
		// Get all samples from total abundance map
		Set<String> allSamples = new HashSet<>(sampleGenomeTotalAbundance.keySet());
		
		// Define marker categories
		String[] markerCategories = {"Anchor", "HighRes", "Growth", "Housekeeping"};
		
		for (String sample : allSamples) {
			Map<String, Double> genomeTotals = sampleGenomeTotalAbundance.get(sample);
			Map<String, Map<String, Double>> genomeCategories = sampleGenomeCategoryAbundance.getOrDefault(sample, new HashMap<>());
			
			// Iterate through all selected genomes (those with total abundance)
			for (String genomeId : genomeTotals.keySet()) {
				if (!selectedGenomes.contains(genomeId)) {
					continue; // Skip if not in selected genomes
				}
				
				double totalGenomeAbundance = genomeTotals.get(genomeId);
				Map<String, Double> categoryAbundances = genomeCategories.get(genomeId);
				
				// If genome has category abundances, write them
				if (categoryAbundances != null && !categoryAbundances.isEmpty()) {
					for (Map.Entry<String, Double> categoryEntry : categoryAbundances.entrySet()) {
						String category = categoryEntry.getKey();
						double abundance = categoryEntry.getValue();
						double ratio = totalGenomeAbundance > 0 ? (abundance / totalGenomeAbundance) : 0.0;
						
						writer.printf("%s\t%s\t%s\t%.6f\t%.6f\t%.6f\n",
							sample, genomeId, category, abundance, totalGenomeAbundance, ratio);
					}
				}
				
				// Also write entries for marker categories that this genome doesn't have
				// (to ensure all marker categories appear for all genomes)
				Set<String> existingCategories = categoryAbundances != null ? categoryAbundances.keySet() : new HashSet<>();
				for (String category : markerCategories) {
					if (!existingCategories.contains(category)) {
						// Genome doesn't have this marker category, write 0.0
						writer.printf("%s\t%s\t%s\t%.6f\t%.6f\t%.6f\n",
							sample, genomeId, category, 0.0, totalGenomeAbundance, 0.0);
					}
				}
			}
			
			// Also include selected genomes that have no total abundance (no peptides at all)
			for (String genomeId : selectedGenomes) {
				if (!genomeTotals.containsKey(genomeId)) {
					// Genome has no peptides, write 0.0 for all categories
					for (String category : markerCategories) {
						writer.printf("%s\t%s\t%s\t%.6f\t%.6f\t%.6f\n",
							sample, genomeId, category, 0.0, 0.0, 0.0);
					}
				}
			}
		}
		
		writer.close();
		System.out.println("Results written to: " + outputFile.getPath());
	}
	
	/**
	 * Calculate marker-4 protein abundances using only distinct peptides (unique to one selected genome)
	 * @param proteinToGenome If provided, uses this mapping instead of extracting genome ID from protein name
	 */
	private void calculateMarkerAbundancesDistinct(Map<String, PeptideData> peptideData, File faaFolder,
			HashSet<String> selectedGenomes, Map<String, Set<String>> peptideToGenomes, 
			File outputFile, AnalysisStatistics stats, Map<String, String> proteinToGenome) throws IOException {
		
		System.out.println("Identifying distinct peptides within selected genomes...");
		
		// First, identify distinct peptides (peptides that appear in only ONE selected genome)
		Set<String> distinctPeptides = new HashSet<>();
		for (Map.Entry<String, Set<String>> entry : peptideToGenomes.entrySet()) {
			String peptide = entry.getKey();
			Set<String> allGenomes = entry.getValue();
			
			// Filter to only selected genomes
			Set<String> selectedGenomesForPeptide = new HashSet<>();
			for (String genomeId : allGenomes) {
				if (selectedGenomes.contains(genomeId)) {
					selectedGenomesForPeptide.add(genomeId);
				}
			}
			
			// If peptide appears in exactly one selected genome, it's distinct
			if (selectedGenomesForPeptide.size() == 1) {
				distinctPeptides.add(peptide);
			}
		}
		
		System.out.println("Found " + distinctPeptides.size() + " distinct peptides within selected genomes");
		
		System.out.println("Loading protein-to-category mapping...");
		Map<String, String> proteinToCategory = loadMarkerProteinCategories(faaFolder);
		
		// Map: sample -> genome -> category -> abundance
		Map<String, Map<String, Map<String, Double>>> sampleGenomeCategoryAbundance = new HashMap<>();
		// Map: sample -> genome -> total abundance (all proteins)
		Map<String, Map<String, Double>> sampleGenomeTotalAbundance = new HashMap<>();
		
		System.out.println("Calculating abundances using distinct peptides only...");
		// Calculate abundances using only distinct peptides
		for (Map.Entry<String, PeptideData> entry : peptideData.entrySet()) {
			String peptide = entry.getKey();
			
			// Skip if not distinct
			if (!distinctPeptides.contains(peptide)) {
				continue;
			}
			
			PeptideData data = entry.getValue();
			
			// Group proteins by genome (should be only one genome for distinct peptides)
			Map<String, Set<String>> genomeProteins = new HashMap<>();
			for (String proteinId : data.proteinIds) {
				String genomeId;
				if (proteinToGenome != null) {
					// Use provided protein-to-genome mapping
					genomeId = proteinToGenome.get(proteinId);
					if (genomeId == null) {
						continue;  // Skip if protein not found in mapping
					}
				} else {
					// Extract genome ID from protein name
					genomeId = extractGenomeId(proteinId);
				}
				
				if (selectedGenomes.contains(genomeId)) {
					genomeProteins.computeIfAbsent(genomeId, k -> new HashSet<>()).add(proteinId);
				}
			}
			
			if (genomeProteins.isEmpty()) {
				continue;
			}
			
			// Track distinct peptides used in abundance calculation
			stats.distinctPeptidesUsedInAbundance++;
			
			// For each sample
			for (Map.Entry<String, Double> intensityEntry : data.intensities.entrySet()) {
				String sample = intensityEntry.getKey();
				double intensity = intensityEntry.getValue();
				
				// For each genome that has this peptide (should be only one for distinct peptides)
				for (Map.Entry<String, Set<String>> genomeEntry : genomeProteins.entrySet()) {
					String genomeId = genomeEntry.getKey();
					Set<String> proteins = genomeEntry.getValue();
					
					// Distribute intensity equally among proteins in this genome
					double intensityPerProtein = intensity / proteins.size();
					
					// Add to total genome abundance (use full intensity, not per-protein)
					sampleGenomeTotalAbundance.computeIfAbsent(sample, k -> new HashMap<>())
						.merge(genomeId, intensity, Double::sum);
					
					// Add to category abundances (distribute per protein)
					for (String proteinId : proteins) {
						String category = proteinToCategory.get(proteinId);
						if (category != null) {
							sampleGenomeCategoryAbundance.computeIfAbsent(sample, k -> new HashMap<>())
								.computeIfAbsent(genomeId, k -> new HashMap<>())
								.merge(category, intensityPerProtein, Double::sum);
						}
					}
				}
			}
		}
		
		// Write results
		System.out.println("Writing distinct peptide results...");
		PrintWriter writer = new PrintWriter(outputFile);
		writer.println("Sample\tGenome\tCategory\tAbundance\tTotal_Genome_Abundance\tRatio");
		
		// Get all samples from total abundance map
		Set<String> allSamples = new HashSet<>(sampleGenomeTotalAbundance.keySet());
		
		// Define marker categories
		String[] markerCategories = {"Anchor", "HighRes", "Growth", "Housekeeping"};
		
		for (String sample : allSamples) {
			Map<String, Double> genomeTotals = sampleGenomeTotalAbundance.get(sample);
			Map<String, Map<String, Double>> genomeCategories = sampleGenomeCategoryAbundance.getOrDefault(sample, new HashMap<>());
			
			// Iterate through all selected genomes (those with total abundance)
			for (String genomeId : genomeTotals.keySet()) {
				if (!selectedGenomes.contains(genomeId)) {
					continue; // Skip if not in selected genomes
				}
				
				double totalGenomeAbundance = genomeTotals.get(genomeId);
				Map<String, Double> categoryAbundances = genomeCategories.get(genomeId);
				
				// If genome has category abundances, write them
				if (categoryAbundances != null && !categoryAbundances.isEmpty()) {
					for (Map.Entry<String, Double> categoryEntry : categoryAbundances.entrySet()) {
						String category = categoryEntry.getKey();
						double abundance = categoryEntry.getValue();
						double ratio = totalGenomeAbundance > 0 ? (abundance / totalGenomeAbundance) : 0.0;
						
						writer.printf("%s\t%s\t%s\t%.6f\t%.6f\t%.6f\n",
							sample, genomeId, category, abundance, totalGenomeAbundance, ratio);
					}
				}
				
				// Also write entries for marker categories that this genome doesn't have
				// (to ensure all marker categories appear for all genomes)
				Set<String> existingCategories = categoryAbundances != null ? categoryAbundances.keySet() : new HashSet<>();
				for (String category : markerCategories) {
					if (!existingCategories.contains(category)) {
						// Genome doesn't have this marker category, write 0.0
						writer.printf("%s\t%s\t%s\t%.6f\t%.6f\t%.6f\n",
							sample, genomeId, category, 0.0, totalGenomeAbundance, 0.0);
					}
				}
			}
			
			// Also include selected genomes that have no total abundance (no peptides at all)
			for (String genomeId : selectedGenomes) {
				if (!genomeTotals.containsKey(genomeId)) {
					// Genome has no peptides, write 0.0 for all categories
					for (String category : markerCategories) {
						writer.printf("%s\t%s\t%s\t%.6f\t%.6f\t%.6f\n",
							sample, genomeId, category, 0.0, 0.0, 0.0);
					}
				}
			}
		}
		
		writer.close();
		System.out.println("Distinct peptide results written to: " + outputFile.getPath());
	}
	
	/**
	 * Main analysis method (overload without predicted folder - uses protein IDs for mapping)
	 */
	public void analyzeDataset(File peptideFile, File faaFolder, File outputFile, String projectName,
			double razorThreshold, double stepThreshold, int maxGenomes) throws IOException {
		analyzeDataset(peptideFile, faaFolder, null, outputFile, projectName, razorThreshold, stepThreshold, maxGenomes);
	}
	
	/**
	 * Main analysis method
	 * @param predictedFolder If provided, maps peptides to genomes by reading from predicted files.
	 *                        If null, uses protein IDs from peptide data (requires genome ID in protein name).
	 */
	public void analyzeDataset(File peptideFile, File faaFolder, File predictedFolder, File outputFile, String projectName,
			double razorThreshold, double stepThreshold, int maxGenomes) throws IOException {

		AnalysisStatistics stats = new AnalysisStatistics();

		System.out.println("=== Dataset Marker-4 Analysis ===");
		System.out.println("Peptide file: " + peptideFile.getPath());
		System.out.println("FAA folder: " + faaFolder.getPath());
		if (predictedFolder != null) {
			System.out.println("Predicted folder: " + predictedFolder.getPath());
		}

		// Step 1: Read peptides and map to genomes
		System.out.println("\nStep 1: Reading peptides and mapping to genomes...");
		// If using predicted files, don't filter by MGYG (proteins will be filtered
		// when mapping from predicted files)
		// If using protein ID extraction, filter by MGYG (requires genome ID in protein
		// name)
		boolean filterMGYGOnly = (predictedFolder == null);
		Map<String, PeptideData> peptideData = readPeptideIntensities(peptideFile, stats, filterMGYGOnly);
		if (filterMGYGOnly) {
			System.out.println("Found " + peptideData.size() + " peptides with MGYG proteins and intensities");
		} else {
			System.out.println("Found " + peptideData.size() + " peptides with valid proteins and intensities");
		}

		Map<String, Set<String>> peptideToGenomes;
		Map<String, String> proteinToGenome = null; // For predicted file mapping
		Map<String, Double> peptideToDetectability = null; // Detectability from predicted files
		if (predictedFolder != null) {
			// Use predicted files to map peptides to genomes
			PredictedMappingResult mappingResult = mapPeptidesToGenomesFromPredictedFiles(predictedFolder, peptideData);
			peptideToGenomes = mappingResult.peptideToGenomes;
			proteinToGenome = mappingResult.proteinToGenome;
			peptideToDetectability = mappingResult.peptideToDetectability;
		} else {
			// Use protein IDs to map peptides to genomes (requires genome ID in protein
			// name)
			peptideToGenomes = mapPeptidesToGenomes(peptideData);
			// Calculate detectability as number of genomes (fallback when predicted files not available)
			peptideToDetectability = new HashMap<>();
			for (Map.Entry<String, Set<String>> entry : peptideToGenomes.entrySet()) {
				peptideToDetectability.put(entry.getKey(), (double) entry.getValue().size());
			}
		}
		stats.peptidesMappedToGenomes = peptideToGenomes.size();
		System.out.println("Mapped " + peptideToGenomes.size() + " peptides to genomes");

		// Step 2: Build genome-to-peptides map
		System.out.println("\nStep 2: Building genome-to-peptides map...");
		HashMap<String, HashSet<String>> genomeToPeptides = new HashMap<>();
		Set<String> allPeptides = new HashSet<>(peptideToGenomes.keySet());
		for (Map.Entry<String, Set<String>> entry : peptideToGenomes.entrySet()) {
			String peptide = entry.getKey();
			for (String genomeId : entry.getValue()) {
				genomeToPeptides.computeIfAbsent(genomeId, k -> new HashSet<>()).add(peptide);
			}
		}
		stats.totalGenomesFound = genomeToPeptides.size();
		System.out.println("Found " + genomeToPeptides.size() + " genomes");

		// Step 3: Apply Occam's Razor
		System.out.println("\nStep 3: Applying Occam's Razor (threshold=" + razorThreshold + ")...");
		File rankingFile = new File(new File(outputFile, "genome_ranks"), projectName + "_genome_ranking.tsv");
		HashSet<String> selectedGenomes = applyOccamsRazor(genomeToPeptides, allPeptides, razorThreshold, stepThreshold,
				maxGenomes, rankingFile);
		stats.genomesSelectedByRazor = selectedGenomes.size();
		stats.peptidesUsedInRazor = allPeptides.size();
		System.out.println("Selected " + selectedGenomes.size() + " genomes");

		// Step 4: Calculate marker-4 abundances (all peptides)
		System.out.println("\nStep 4: Calculating marker-4 protein abundances (all peptides)...");
		File abundanceOutputFile = new File(new File(outputFile, "Figure1D"), projectName + "_marker_abundances.tsv");
		calculateMarkerAbundances(peptideData, faaFolder, selectedGenomes, abundanceOutputFile, stats, proteinToGenome);

		// Step 4b: Calculate marker-4 abundances (distinct peptides only)
		System.out.println("\nStep 4b: Calculating marker-4 protein abundances (distinct peptides only)...");
		File distinctOutputFile = new File(new File(outputFile, "Figure1E"),
				projectName + "_marker_abundances_unique.tsv");
		calculateMarkerAbundancesDistinct(peptideData, faaFolder, selectedGenomes, peptideToGenomes, distinctOutputFile,
				stats, proteinToGenome);

		// Step 5: Write peptide detectability statistics
		System.out.println("\nStep 5: Writing peptide detectability statistics...");
		File detectabilityFile = new File(new File(outputFile, "Figure1F"), projectName + "_peptide_detectability.tsv");
		writePeptideDetectability(peptideData, peptideToDetectability, faaFolder, detectabilityFile);
		System.out.println("Peptide detectability written to: " + detectabilityFile.getPath());

		// Step 6: Write general statistics
		System.out.println("\nStep 6: Writing general statistics...");
		File statsFile = new File(new File(outputFile, "statistics"), projectName + "_statistics.txt");
		stats.writeToFile(statsFile, peptideToGenomes);
		System.out.println("Statistics written to: " + statsFile.getPath());

		System.out.println("\n=== Analysis Complete ===");
	}
	
	/**
	 * Analyze dataset using protein info file (from extractProteinInfoFromFasta) instead of predicted folder
	 * This version gets marker protein information from the protein info file and doesn't calculate detectability
	 * Uses identifier (first part before underscore in reference, e.g., "ATN" from "ATN_Q7CUX4") as genome identifier
	 * 
	 * @param peptideFile Input peptide-intensity TSV file
	 * @param proteinInfoFile Protein info TSV file (Reference, Taxa, Marker_Category) from extractProteinInfoFromFasta
	 * @param outputFile Output directory
	 * @param projectName Project name for output files
	 * @param razorThreshold Threshold for Occam's Razor (e.g., 0.95)
	 * @param stepThreshold Step threshold for Occam's Razor (e.g., 0.00001)
	 * @param maxGenomes Maximum genomes to select
	 */
	public void analyzeDatasetFromProteinInfo(File peptideFile, File proteinInfoFile, File outputFile, String projectName,
			double razorThreshold, double stepThreshold, int maxGenomes) throws IOException {
		
		AnalysisStatistics stats = new AnalysisStatistics();
		
		System.out.println("=== Dataset Marker-4 Analysis (from Protein Info) ===");
		System.out.println("Peptide file: " + peptideFile.getPath());
		System.out.println("Protein info file: " + proteinInfoFile.getPath());
		
		// Step 1: Load protein info file
		System.out.println("\nStep 1: Loading protein info file...");
		Map<String, String> referenceToTaxa = new HashMap<>();
		Map<String, String> referenceToCategory = new HashMap<>();
		Map<String, String> proteinToGenome = new HashMap<>();
		
		try (BufferedReader reader = new BufferedReader(new FileReader(proteinInfoFile))) {
			String header = reader.readLine(); // Skip header
			if (header == null || !header.contains("Reference")) {
				throw new IOException("Invalid protein info file format. Expected header: Reference\tTaxa\tMarker_Category");
			}
			
			String line;
			int lineNum = 1;
			while ((line = reader.readLine()) != null) {
				lineNum++;
				String[] parts = line.split("\t");
				if (parts.length < 3) {
					System.err.println("Warning: Skipping line " + lineNum + " (invalid format): " + line);
					continue;
				}
				
				String reference = parts[0].trim();
				String taxa = parts[1].trim();
				String category = parts[2].trim();
				
				// Extract identifier (first part before underscore) as genome identifier
				String genome = null;
				int underscoreIndex = reference.indexOf('_');
				if (underscoreIndex > 0) {
					genome = reference.substring(0, underscoreIndex);
				} else {
					// No underscore, use entire reference as identifier
					genome = reference;
				}
				
				if (!reference.isEmpty() && !taxa.isEmpty() && genome != null && !genome.isEmpty()) {
					referenceToTaxa.put(reference, taxa);
					proteinToGenome.put(reference, genome); // Use identifier (e.g., "ATN") as genome identifier
					if (!category.equals("None") && !category.isEmpty()) {
						referenceToCategory.put(reference, category);
					}
				}
			}
		}
		
		System.out.println("Loaded " + referenceToTaxa.size() + " protein references");
		System.out.println("Found " + referenceToCategory.size() + " marker proteins");
		
		// Step 2: Read peptides and map to genomes using protein info
		System.out.println("\nStep 2: Reading peptides and mapping to genomes...");
		Map<String, PeptideData> peptideData = readPeptideIntensities(peptideFile, stats, false);
		System.out.println("Found " + peptideData.size() + " peptides with valid proteins and intensities");
		
		// Map peptides to genomes using protein references (use identifier as genome)
		Map<String, Set<String>> peptideToGenomes = new HashMap<>();
		for (Map.Entry<String, PeptideData> entry : peptideData.entrySet()) {
			String peptide = entry.getKey();
			PeptideData data = entry.getValue();

			Set<String> genomes = new HashSet<>();
			for (String proteinId : data.proteinIds) {
				// Use proteinToGenome map to get identifier (genome)
				String genome = proteinToGenome.get(proteinId);
				if (genome != null && !genome.isEmpty()) {
					genomes.add(genome);
				} else {
					// Fallback: extract identifier directly from proteinId if not in map
					int underscoreIndex = proteinId.indexOf('_');
					if (underscoreIndex > 0) {
						genome = proteinId.substring(0, underscoreIndex);
						if (!genome.isEmpty()) {
							genomes.add(genome);
						}
					} else {
						// No underscore, use entire proteinId as identifier
						genomes.add(proteinId);
					}
				}
			}

			if (!genomes.isEmpty()) {
				peptideToGenomes.put(peptide, genomes);
			}
		}
		
		stats.peptidesMappedToGenomes = peptideToGenomes.size();
		System.out.println("Mapped " + peptideToGenomes.size() + " peptides to genomes");
		
		// Step 3: Build genome-to-peptides map
		System.out.println("\nStep 3: Building genome-to-peptides map...");
		HashMap<String, HashSet<String>> genomeToPeptides = new HashMap<>();
		Set<String> allPeptides = new HashSet<>(peptideToGenomes.keySet());
		for (Map.Entry<String, Set<String>> entry : peptideToGenomes.entrySet()) {
			String peptide = entry.getKey();
			for (String genomeId : entry.getValue()) {
				genomeToPeptides.computeIfAbsent(genomeId, k -> new HashSet<>()).add(peptide);
			}
		}
		stats.totalGenomesFound = genomeToPeptides.size();
		System.out.println("Found " + genomeToPeptides.size() + " genomes");
		
		// Step 4: Apply Occam's Razor
		System.out.println("\nStep 4: Applying Occam's Razor (threshold=" + razorThreshold + ")...");
		File rankingFile = new File(new File(outputFile, "genome_ranks"), projectName + "_genome_ranking.tsv");
		HashSet<String> selectedGenomes = applyOccamsRazor(genomeToPeptides, allPeptides, razorThreshold, stepThreshold,
				maxGenomes, rankingFile);
		stats.genomesSelectedByRazor = selectedGenomes.size();
		stats.peptidesUsedInRazor = allPeptides.size();
		System.out.println("Selected " + selectedGenomes.size() + " genomes");
		
		// Step 5: Calculate marker-4 abundances (all peptides)
		System.out.println("\nStep 5: Calculating marker-4 protein abundances (all peptides)...");
		File abundanceOutputFile = new File(new File(outputFile, "Figure1D"), projectName + "_marker_abundances.tsv");
		calculateMarkerAbundancesFromProteinInfo(peptideData, referenceToCategory, proteinToGenome, 
				selectedGenomes, abundanceOutputFile, stats, peptideToGenomes, false);
		
		// Step 5b: Calculate marker-4 abundances (distinct peptides only)
		System.out.println("\nStep 5b: Calculating marker-4 protein abundances (distinct peptides only)...");
		File distinctOutputFile = new File(new File(outputFile, "Figure1E"), 
				projectName + "_marker_abundances_unique.tsv");
		calculateMarkerAbundancesFromProteinInfo(peptideData, referenceToCategory, proteinToGenome, 
				selectedGenomes, distinctOutputFile, stats, peptideToGenomes, true);
		
		// Step 5c: Calculate genome abundance using median of peptides (10 methods)
		System.out.println("\nStep 5c: Calculating genome abundance using median of peptides...");
		File medianOutputFile = new File(new File(outputFile, "Figure2"), 
				projectName + "_genome_abundance_median.tsv");
		calculateGenomeAbundanceMedian(peptideData, referenceToCategory, proteinToGenome, 
				selectedGenomes, medianOutputFile, peptideToGenomes);
		
		// Step 6: Write general statistics
		System.out.println("\nStep 6: Writing general statistics...");
		File statsFile = new File(new File(outputFile, "statistics"), projectName + "_statistics.txt");
		stats.writeToFile(statsFile, peptideToGenomes);
		System.out.println("Statistics written to: " + statsFile.getPath());
		
		System.out.println("\n=== Analysis Complete ===");
	}
	
	/**
	 * Calculate marker-4 protein abundances using protein info file instead of FAA folder
	 * Uses identifier (first part before underscore) as genome identifier
	 * @param distinctOnly If true, only use distinct peptides (unique to one selected genome)
	 */
	private void calculateMarkerAbundancesFromProteinInfo(Map<String, PeptideData> peptideData,
			Map<String, String> referenceToCategory, Map<String, String> proteinToGenome,
			HashSet<String> selectedGenomes, File outputFile, AnalysisStatistics stats,
			Map<String, Set<String>> peptideToGenomes, boolean distinctOnly) throws IOException {
		
		System.out.println("Selected genome count: " + selectedGenomes.size());
		
		// If distinct only, identify distinct peptides first
		Set<String> distinctPeptides = null;
		if (distinctOnly) {
			System.out.println("Identifying distinct peptides within selected genomes...");
			distinctPeptides = new HashSet<>();
			for (Map.Entry<String, Set<String>> entry : peptideToGenomes.entrySet()) {
				String peptide = entry.getKey();
				Set<String> allGenomes = entry.getValue();
				
				// Filter to only selected genomes
				Set<String> selectedGenomesForPeptide = new HashSet<>();
				for (String genomeId : allGenomes) {
					if (selectedGenomes.contains(genomeId)) {
						selectedGenomesForPeptide.add(genomeId);
					}
				}
				
				// If peptide appears in exactly one selected genome, it's distinct
				if (selectedGenomesForPeptide.size() == 1) {
					distinctPeptides.add(peptide);
				}
			}
			System.out.println("Found " + distinctPeptides.size() + " distinct peptides within selected genomes");
		}
		
		// Map: sample -> genome -> category -> abundance
		Map<String, Map<String, Map<String, Double>>> sampleGenomeCategoryAbundance = new HashMap<>();
		// Map: sample -> genome -> total abundance (all proteins)
		Map<String, Map<String, Double>> sampleGenomeTotalAbundance = new HashMap<>();
		
		System.out.println("Calculating abundances" + (distinctOnly ? " using distinct peptides only" : "") + "...");
		
		// Calculate abundances
		for (Map.Entry<String, PeptideData> entry : peptideData.entrySet()) {
			String peptide = entry.getKey();
			
			// Skip if distinct only and not distinct
			if (distinctOnly && (distinctPeptides == null || !distinctPeptides.contains(peptide))) {
				continue;
			}
			
			PeptideData data = entry.getValue();
			
			// Group proteins by genome (using identifier)
			Map<String, Set<String>> genomeProteins = new HashMap<>();
			for (String proteinId : data.proteinIds) {
				String genomeId = proteinToGenome.get(proteinId);
				if (genomeId == null || !selectedGenomes.contains(genomeId)) {
					// Fallback: extract identifier directly if not in map
					int underscoreIndex = proteinId.indexOf('_');
					if (underscoreIndex > 0) {
						genomeId = proteinId.substring(0, underscoreIndex);
					} else {
						genomeId = proteinId;
					}
					if (!selectedGenomes.contains(genomeId)) {
						continue;
					}
				}
				// Skip excluded genomes
				if (EXCLUDED_GENOMES.contains(genomeId)) {
					continue;
				}
				genomeProteins.computeIfAbsent(genomeId, k -> new HashSet<>()).add(proteinId);
			}
			
			if (genomeProteins.isEmpty()) {
				continue;
			}
			
			// Track peptides used in abundance calculation
			if (distinctOnly) {
				stats.distinctPeptidesUsedInAbundance++;
			} else {
				stats.peptidesUsedInAbundance++;
			}
			
			// For each sample
			for (Map.Entry<String, Double> intensityEntry : data.intensities.entrySet()) {
				String sample = intensityEntry.getKey();
				double intensity = intensityEntry.getValue();
				
				// For each genome that has this peptide
				for (Map.Entry<String, Set<String>> genomeEntry : genomeProteins.entrySet()) {
					String genomeId = genomeEntry.getKey();
					Set<String> proteins = genomeEntry.getValue();
					
					// Distribute intensity equally among proteins in this genome
					// (razor principle: if peptide maps to multiple proteins in same genome, split equally)
					double intensityPerProtein = intensity / proteins.size();
					
					// Add to total genome abundance (use full intensity, not per-protein)
					sampleGenomeTotalAbundance.computeIfAbsent(sample, k -> new HashMap<>())
						.merge(genomeId, intensity, Double::sum);
					
					// Add to category abundances (distribute per protein)
					for (String proteinId : proteins) {
						String category = referenceToCategory.get(proteinId);
						if (category != null) {
							sampleGenomeCategoryAbundance.computeIfAbsent(sample, k -> new HashMap<>())
								.computeIfAbsent(genomeId, k -> new HashMap<>())
								.merge(category, intensityPerProtein, Double::sum);
						}
					}
				}
			}
		}
		
		// Write results
		System.out.println("Writing results...");
		PrintWriter writer = new PrintWriter(outputFile);
		writer.println("Sample\tGenome\tCategory\tAbundance\tTotal_Genome_Abundance\tRatio");
		
		// Get all samples from total abundance map
		Set<String> allSamples = new HashSet<>(sampleGenomeTotalAbundance.keySet());
		
		// Define marker categories
		String[] markerCategories = {"Anchor", "HighRes", "Growth", "Housekeeping"};
		
		for (String sample : allSamples) {
			Map<String, Double> genomeTotals = sampleGenomeTotalAbundance.get(sample);
			Map<String, Map<String, Double>> genomeCategories = sampleGenomeCategoryAbundance.getOrDefault(sample, new HashMap<>());
			
			// Iterate through all selected genomes (those with total abundance)
			for (String genomeId : genomeTotals.keySet()) {
				if (!selectedGenomes.contains(genomeId)) {
					continue; // Skip if not in selected genomes
				}
				// Skip excluded genomes (no function annotation)
				if (EXCLUDED_GENOMES.contains(genomeId)) {
					continue;
				}
				
				double totalGenomeAbundance = genomeTotals.get(genomeId);
				Map<String, Double> categoryAbundances = genomeCategories.get(genomeId);
				
				// If genome has category abundances, write them
				if (categoryAbundances != null && !categoryAbundances.isEmpty()) {
					for (Map.Entry<String, Double> categoryEntry : categoryAbundances.entrySet()) {
						String category = categoryEntry.getKey();
						double abundance = categoryEntry.getValue();
						double ratio = totalGenomeAbundance > 0 ? (abundance / totalGenomeAbundance) : 0.0;
						
						writer.printf("%s\t%s\t%s\t%.6f\t%.6f\t%.6f\n",
							sample, genomeId, category, abundance, totalGenomeAbundance, ratio);
					}
				}
				
				// Also write entries for marker categories that this genome doesn't have
				// (to ensure all marker categories appear for all genomes)
				Set<String> existingCategories = categoryAbundances != null ? categoryAbundances.keySet() : new HashSet<>();
				for (String category : markerCategories) {
					if (!existingCategories.contains(category)) {
						// Genome doesn't have this marker category, write 0.0
						writer.printf("%s\t%s\t%s\t%.6f\t%.6f\t%.6f\n",
							sample, genomeId, category, 0.0, totalGenomeAbundance, 0.0);
					}
				}
			}
			
			// Also include selected genomes that have no total abundance (no peptides at all)
			for (String genomeId : selectedGenomes) {
				// Skip excluded genomes (no function annotation)
				if (EXCLUDED_GENOMES.contains(genomeId)) {
					continue;
				}
				if (!genomeTotals.containsKey(genomeId)) {
					// Genome has no peptides, write 0.0 for all categories
					for (String category : markerCategories) {
						writer.printf("%s\t%s\t%s\t%.6f\t%.6f\t%.6f\n",
							sample, genomeId, category, 0.0, 0.0, 0.0);
					}
				}
			}
		}
		
		writer.close();
		String outputType = distinctOnly ? "distinct peptide" : "";
		System.out.println("Marker abundance " + outputType + " results written to: " + outputFile.getPath());
	}
	
	/**
	 * Calculate genome abundance using median of peptides (10 different methods)
	 * Methods:
	 * 1. Median of Anchor peptides (all)
	 * 2. Median of HighRes peptides (all)
	 * 3. Median of Growth peptides (all)
	 * 4. Median of Housekeeping peptides (all)
	 * 5. Median of all peptides
	 * 6. Median of Anchor peptides (distinct only)
	 * 7. Median of HighRes peptides (distinct only)
	 * 8. Median of Growth peptides (distinct only)
	 * 9. Median of Housekeeping peptides (distinct only)
	 * 10. Median of all distinct peptides
	 */
	private void calculateGenomeAbundanceMedian(Map<String, PeptideData> peptideData,
			Map<String, String> referenceToCategory, Map<String, String> proteinToGenome,
			HashSet<String> selectedGenomes, File outputFile, Map<String, Set<String>> peptideToGenomes) throws IOException {
		
		System.out.println("Identifying distinct peptides within selected genomes...");
		
		// Identify distinct peptides (unique to one selected genome)
		Set<String> distinctPeptides = new HashSet<>();
		for (Map.Entry<String, Set<String>> entry : peptideToGenomes.entrySet()) {
			String peptide = entry.getKey();
			Set<String> allGenomes = entry.getValue();
			
			// Filter to only selected genomes
			Set<String> selectedGenomesForPeptide = new HashSet<>();
			for (String genomeId : allGenomes) {
				if (selectedGenomes.contains(genomeId) && !EXCLUDED_GENOMES.contains(genomeId)) {
					selectedGenomesForPeptide.add(genomeId);
				}
			}
			
			// If peptide appears in exactly one selected genome, it's distinct
			if (selectedGenomesForPeptide.size() == 1) {
				distinctPeptides.add(peptide);
			}
		}
		
		System.out.println("Found " + distinctPeptides.size() + " distinct peptides within selected genomes");
		
		// Map peptides to their categories
		Map<String, Set<String>> peptideToCategories = new HashMap<>();
		for (Map.Entry<String, PeptideData> entry : peptideData.entrySet()) {
			String peptide = entry.getKey();
			PeptideData data = entry.getValue();
			Set<String> categories = new HashSet<>();
			
			for (String proteinId : data.proteinIds) {
				String category = referenceToCategory.get(proteinId);
				if (category != null && !category.isEmpty()) {
					categories.add(category);
				}
			}
			
			if (!categories.isEmpty()) {
				peptideToCategories.put(peptide, categories);
			}
		}
		
		// Define marker categories
		String[] markerCategories = {"Anchor", "HighRes", "Growth", "Housekeeping"};
		
		// Structure: sample -> genome -> method -> median
		// Methods: "Anchor_All", "HighRes_All", "Growth_All", "Housekeeping_All", "All_All",
		//          "Anchor_Distinct", "HighRes_Distinct", "Growth_Distinct", "Housekeeping_Distinct", "All_Distinct"
		Map<String, Map<String, Map<String, Double>>> results = new HashMap<>();
		
		// Get all samples from peptide data
		Set<String> allSamples = new HashSet<>();
		for (PeptideData data : peptideData.values()) {
			allSamples.addAll(data.intensities.keySet());
		}
		
		System.out.println("Calculating medians for " + allSamples.size() + " samples...");
		
		// For each sample and genome
		for (String sample : allSamples) {
			results.put(sample, new HashMap<>());
			
			for (String genomeId : selectedGenomes) {
				// Skip excluded genomes
				if (EXCLUDED_GENOMES.contains(genomeId)) {
					continue;
				}
				
				results.get(sample).put(genomeId, new HashMap<>());
				
				// Collect intensities for each method
				Map<String, ArrayList<Double>> methodIntensities = new HashMap<>();
				
				for (String method : new String[]{"Anchor_All", "HighRes_All", "Growth_All", "Housekeeping_All", "All_All",
						"Anchor_Distinct", "HighRes_Distinct", "Growth_Distinct", "Housekeeping_Distinct", "All_Distinct"}) {
					methodIntensities.put(method, new ArrayList<>());
				}
				
				// Process each peptide
				for (Map.Entry<String, PeptideData> entry : peptideData.entrySet()) {
					String peptide = entry.getKey();
					PeptideData data = entry.getValue();
					
					// Check if peptide maps to this genome
					Set<String> peptideGenomes = peptideToGenomes.get(peptide);
					if (peptideGenomes == null || !peptideGenomes.contains(genomeId)) {
						continue;
					}
					
					// Get intensity for this sample
					Double intensity = data.intensities.get(sample);
					if (intensity == null || intensity <= 0) {
						continue;
					}
					
					// Check if peptide is distinct for this genome
					// A peptide is distinct if it's in the distinctPeptides set (appears in exactly one selected genome)
					// Since we already verified peptideGenomes.contains(genomeId), if it's in distinctPeptides, 
					// it must be distinct for this genome
					boolean isDistinct = distinctPeptides.contains(peptide);
					
					// Get categories for this peptide
					Set<String> categories = peptideToCategories.get(peptide);
					boolean hasCategory = categories != null && !categories.isEmpty();
					
					// Add to "All_All" (all peptides)
					methodIntensities.get("All_All").add(intensity);
					
					// Add to "All_Distinct" if distinct
					if (isDistinct) {
						methodIntensities.get("All_Distinct").add(intensity);
					}
					
					// Add to category-specific methods
					if (hasCategory) {
						for (String category : categories) {
							// Add to category_All
							methodIntensities.get(category + "_All").add(intensity);
							
							// Add to category_Distinct if distinct
							if (isDistinct) {
								methodIntensities.get(category + "_Distinct").add(intensity);
							}
						}
					}
				}
				
				// Calculate median for each method
				for (Map.Entry<String, ArrayList<Double>> methodEntry : methodIntensities.entrySet()) {
					String method = methodEntry.getKey();
					ArrayList<Double> intensities = methodEntry.getValue();
					
					if (!intensities.isEmpty()) {
						double median = calculateMedian(intensities);
						results.get(sample).get(genomeId).put(method, median);
					}
				}
			}
		}
		
		// Write results
		System.out.println("Writing median abundance results...");
		PrintWriter writer = new PrintWriter(outputFile);
		
		// Header: Sample, Genome, Anchor_All, HighRes_All, Growth_All, Housekeeping_All, All_All,
		//         Anchor_Distinct, HighRes_Distinct, Growth_Distinct, Housekeeping_Distinct, All_Distinct
		writer.println("Sample\tGenome\tAnchor_All\tHighRes_All\tGrowth_All\tHousekeeping_All\tAll_All\t" +
				"Anchor_Distinct\tHighRes_Distinct\tGrowth_Distinct\tHousekeeping_Distinct\tAll_Distinct");
		
		for (String sample : allSamples) {
			Map<String, Map<String, Double>> genomeResults = results.get(sample);
			
			for (String genomeId : selectedGenomes) {
				// Skip excluded genomes
				if (EXCLUDED_GENOMES.contains(genomeId)) {
					continue;
				}
				
				Map<String, Double> methodResults = genomeResults.get(genomeId);
				if (methodResults == null) {
					methodResults = new HashMap<>();
				}
				
				writer.printf("%s\t%s", sample, genomeId);
				
				// Write all 10 methods
				for (String method : new String[]{"Anchor_All", "HighRes_All", "Growth_All", "Housekeeping_All", "All_All",
						"Anchor_Distinct", "HighRes_Distinct", "Growth_Distinct", "Housekeeping_Distinct", "All_Distinct"}) {
					Double median = methodResults.get(method);
					if (median != null) {
						writer.printf("\t%.6f", median);
					} else {
						writer.printf("\tNA");
					}
				}
				
				writer.println();
			}
		}
		
		writer.close();
		System.out.println("Median abundance results written to: " + outputFile.getPath());
	}
	
	/**
	 * Calculate median of a list of doubles
	 */
	private double calculateMedian(ArrayList<Double> values) {
		if (values.isEmpty()) {
			return 0.0;
		}
		
		// Sort the values
		Collections.sort(values);
		
		int size = values.size();
		if (size % 2 == 0) {
			// Even number of values - average of two middle values
			double mid1 = values.get(size / 2 - 1);
			double mid2 = values.get(size / 2);
			return (mid1 + mid2) / 2.0;
		} else {
			// Odd number of values - middle value
			return values.get(size / 2);
		}
	}
	
	/**
	 * Calculate relative abundance of taxa from peptide intensities using different methods
	 * and compare with theoretical values from paper_data.txt
	 */
	public void calculateTaxaAbundanceFromPeptides(File peptideFile, File paperDataFile, File outputFile) throws IOException {
		if (!peptideFile.exists() || !peptideFile.isFile()) {
			throw new IOException("Peptide file does not exist: " + peptideFile.getAbsolutePath());
		}
		if (!paperDataFile.exists() || !paperDataFile.isFile()) {
			throw new IOException("Paper data file does not exist: " + paperDataFile.getAbsolutePath());
		}
		
		System.out.println("=== Taxa Abundance Analysis from Peptide Intensities ===");
		System.out.println("Reading peptide file: " + peptideFile.getPath());
		
		// Step 1: Read peptide intensities and map to taxa
		Map<String, Map<String, ArrayList<Double>>> taxaSampleIntensities = new HashMap<>(); // taxa -> sample -> list of intensities
		Map<String, Set<String>> peptideToTaxa = new HashMap<>(); // For distinct peptide calculation
		Map<String, Integer> sampleColIndices = new HashMap<>(); // sample name -> column index (used later)
		
		try (BufferedReader reader = new BufferedReader(new FileReader(peptideFile))) {
			String headerLine = reader.readLine();
			if (headerLine == null) {
				throw new IOException("Empty peptide file");
			}
			
			String[] headers = headerLine.split("\t");
			int peptideIdx = -1;
			int proteinsIdx = -1;
			
			// Find column indices
			for (int i = 0; i < headers.length; i++) {
				String h = headers[i].trim();
				if (h.equalsIgnoreCase("peptide") || h.equalsIgnoreCase("sequence")) {
					peptideIdx = i;
				} else if (h.equalsIgnoreCase("proteins")) {
					proteinsIdx = i;
				} else if (h.contains("Run1_") && h.contains(".mzML")) {
					// Extract sample name: Run1_C1_2000ng.mzML -> C1
					String sampleName = extractSampleName(h);
					if (sampleName != null) {
						sampleColIndices.put(sampleName, i);
					}
				}
			}
			
			if (peptideIdx < 0 || proteinsIdx < 0) {
				throw new IOException("Could not find peptide or proteins column");
			}
			
			if (sampleColIndices.isEmpty()) {
				throw new IOException("Could not find any sample intensity columns");
			}
			
			System.out.println("Found " + sampleColIndices.size() + " sample columns: " + sampleColIndices.keySet());
			
			// Read data
			String line;
			int lineCount = 0;
			while ((line = reader.readLine()) != null) {
				lineCount++;
				if (lineCount % 10000 == 0) {
					System.out.println("Processed " + lineCount + " peptides...");
				}
				
				String[] parts = line.split("\t");
				if (parts.length <= Math.max(peptideIdx, proteinsIdx)) {
					continue;
				}
				
				String peptide = parts[peptideIdx].trim();
				String proteins = parts[proteinsIdx].trim();
				
				if (peptide.isEmpty() || proteins.isEmpty()) {
					continue;
				}
				
				// Extract taxa identifiers from protein names
				Set<String> taxaSet = new HashSet<>();
				for (String protein : proteins.split(";")) {
					protein = protein.trim();
					if (protein.isEmpty()) {
						continue;
					}
					
					// Skip decoys and host proteins
					if (protein.startsWith("REV") || protein.startsWith("sp") || protein.startsWith("CRAP")) {
						continue;
					}
					
					// Extract taxa identifier (part before first underscore)
					int underscoreIdx = protein.indexOf('_');
					if (underscoreIdx > 0) {
						String taxa = protein.substring(0, underscoreIdx);
						if (!taxa.isEmpty()) {
							taxaSet.add(taxa);
						}
					}
				}
				
				if (taxaSet.isEmpty()) {
					continue;
				}
				
				// Track peptide-to-taxa mapping for distinct peptide calculation
				peptideToTaxa.put(peptide, taxaSet);
				
				// Extract intensities for each sample
				for (Map.Entry<String, Integer> sampleEntry : sampleColIndices.entrySet()) {
					String sample = sampleEntry.getKey();
					int colIdx = sampleEntry.getValue();
					
					if (colIdx < parts.length) {
						try {
							double intensity = Double.parseDouble(parts[colIdx].trim());
							if (intensity > 0) {
								// Add intensity to all taxa this peptide maps to
								for (String taxa : taxaSet) {
									taxaSampleIntensities.computeIfAbsent(taxa, k -> new HashMap<>())
										.computeIfAbsent(sample, k -> new ArrayList<>())
										.add(intensity);
								}
							}
						} catch (NumberFormatException e) {
							// Skip invalid numbers
						}
					}
				}
			}
			
			System.out.println("Found " + taxaSampleIntensities.size() + " taxa with peptide data");
		}
		
		// Step 2: Calculate relative abundances using different methods
		System.out.println("\nCalculating relative abundances using different methods...");
		
		// Methods: Sum, Median, Mean, Max, Sum_Distinct
		Map<String, Map<String, Map<String, Double>>> abundances = new HashMap<>(); // method -> sample -> taxa -> abundance
		Map<String, Map<String, Map<String, Integer>>> peptideCounts = new HashMap<>(); // method -> sample -> taxa -> peptide count
		
		for (String method : Arrays.asList("Sum", "Median", "Mean", "Max", "Sum_Distinct")) {
			abundances.put(method, new HashMap<>());
			peptideCounts.put(method, new HashMap<>());
		}
		
		// Identify distinct peptides (unique to one taxa)
		Set<String> distinctPeptides = new HashSet<>();
		for (Map.Entry<String, Set<String>> entry : peptideToTaxa.entrySet()) {
			if (entry.getValue().size() == 1) {
				distinctPeptides.add(entry.getKey());
			}
		}
		System.out.println("Found " + distinctPeptides.size() + " distinct peptides (unique to one taxa)");
		
		// Collect distinct peptide intensities (already collected in first pass, just filter)
		Map<String, Map<String, ArrayList<Double>>> taxaSampleDistinctIntensities = new HashMap<>();
		for (Map.Entry<String, Set<String>> entry : peptideToTaxa.entrySet()) {
			String peptide = entry.getKey();
			Set<String> taxaSet = entry.getValue();
			
			// Only process distinct peptides (unique to one taxa)
			if (taxaSet.size() == 1) {
				String taxa = taxaSet.iterator().next();
				Map<String, ArrayList<Double>> sampleIntensities = taxaSampleIntensities.get(taxa);
				if (sampleIntensities != null) {
					// We need to track which intensities belong to distinct peptides
					// For now, we'll re-read to get distinct peptide intensities separately
					// This is simpler than tracking during first pass
				}
			}
		}
		
		// Re-read file to collect distinct peptide intensities
		if (!distinctPeptides.isEmpty()) {
			try (BufferedReader reader = new BufferedReader(new FileReader(peptideFile))) {
				String headerLine = reader.readLine(); // Skip header
				if (headerLine == null) {
					throw new IOException("Empty peptide file");
				}
				
				String[] headers = headerLine.split("\t");
				int peptideIdx = -1;
				int proteinsIdx = -1;
				
				// Find column indices
				for (int i = 0; i < headers.length; i++) {
					String h = headers[i].trim();
					if (h.equalsIgnoreCase("peptide") || h.equalsIgnoreCase("sequence")) {
						peptideIdx = i;
					} else if (h.equalsIgnoreCase("proteins")) {
						proteinsIdx = i;
					}
				}
				
				if (peptideIdx < 0 || proteinsIdx < 0) {
					throw new IOException("Could not find peptide or proteins column");
				}
				
				String line;
				while ((line = reader.readLine()) != null) {
					String[] parts = line.split("\t");
					if (parts.length <= Math.max(peptideIdx, proteinsIdx)) {
						continue;
					}
					
					String peptide = parts[peptideIdx].trim();
					if (!distinctPeptides.contains(peptide)) {
						continue;
					}
					
					Set<String> taxaSet = peptideToTaxa.get(peptide);
					if (taxaSet == null || taxaSet.size() != 1) {
						continue;
					}
					
					String taxa = taxaSet.iterator().next();
					
					// Extract intensities for each sample
					for (Map.Entry<String, Integer> sampleEntry : sampleColIndices.entrySet()) {
						String sample = sampleEntry.getKey();
						int colIdx = sampleEntry.getValue();
						
						if (colIdx < parts.length) {
							try {
								double intensity = Double.parseDouble(parts[colIdx].trim());
								if (intensity > 0) {
									taxaSampleDistinctIntensities.computeIfAbsent(taxa, k -> new HashMap<>())
										.computeIfAbsent(sample, k -> new ArrayList<>())
										.add(intensity);
								}
							} catch (NumberFormatException e) {
								// Skip
							}
						}
					}
				}
			}
		}
		
		// Calculate abundances for each method and sample
		for (String sample : sampleColIndices.keySet()) {
			// Calculate raw abundances for each method
			Map<String, Double> rawAbundances = new HashMap<>();
			
			for (String taxa : taxaSampleIntensities.keySet()) {
				Map<String, ArrayList<Double>> sampleIntensities = taxaSampleIntensities.get(taxa);
				ArrayList<Double> intensities = sampleIntensities.get(sample);
				
				if (intensities == null || intensities.isEmpty()) {
					continue;
				}
				
				int peptideCount = intensities.size();
				
				// Sum
				double sum = intensities.stream().mapToDouble(Double::doubleValue).sum();
				rawAbundances.put(taxa, sum);
				abundances.get("Sum").computeIfAbsent(sample, k -> new HashMap<>()).put(taxa, sum);
				peptideCounts.get("Sum").computeIfAbsent(sample, k -> new HashMap<>()).put(taxa, peptideCount);
				
				// Median
				double median = calculateMedian(intensities);
				abundances.get("Median").computeIfAbsent(sample, k -> new HashMap<>()).put(taxa, median);
				peptideCounts.get("Median").computeIfAbsent(sample, k -> new HashMap<>()).put(taxa, peptideCount);
				
				// Mean
				double mean = sum / intensities.size();
				abundances.get("Mean").computeIfAbsent(sample, k -> new HashMap<>()).put(taxa, mean);
				peptideCounts.get("Mean").computeIfAbsent(sample, k -> new HashMap<>()).put(taxa, peptideCount);
				
				// Max
				double max = Collections.max(intensities);
				abundances.get("Max").computeIfAbsent(sample, k -> new HashMap<>()).put(taxa, max);
				peptideCounts.get("Max").computeIfAbsent(sample, k -> new HashMap<>()).put(taxa, peptideCount);
			}
			
			// Sum_Distinct
			for (String taxa : taxaSampleDistinctIntensities.keySet()) {
				Map<String, ArrayList<Double>> sampleIntensities = taxaSampleDistinctIntensities.get(taxa);
				ArrayList<Double> intensities = sampleIntensities.get(sample);
				
				if (intensities != null && !intensities.isEmpty()) {
					int distinctPeptideCount = intensities.size();
					double sumDistinct = intensities.stream().mapToDouble(Double::doubleValue).sum();
					abundances.get("Sum_Distinct").computeIfAbsent(sample, k -> new HashMap<>()).put(taxa, sumDistinct);
					peptideCounts.get("Sum_Distinct").computeIfAbsent(sample, k -> new HashMap<>()).put(taxa, distinctPeptideCount);
				}
			}
			
			// Normalize to percentages for each method
			for (String method : Arrays.asList("Sum", "Median", "Mean", "Max", "Sum_Distinct")) {
				Map<String, Double> methodAbundances = abundances.get(method).get(sample);
				if (methodAbundances == null || methodAbundances.isEmpty()) {
					continue;
				}
				
				double total = methodAbundances.values().stream().mapToDouble(Double::doubleValue).sum();
				if (total > 0) {
					for (String taxa : methodAbundances.keySet()) {
						double normalized = (methodAbundances.get(taxa) / total) * 100.0;
						methodAbundances.put(taxa, normalized);
					}
				}
			}
		}
		
		// Step 3: Read theoretical values from paper_data.txt
		System.out.println("Reading theoretical values from paper_data.txt...");
		Map<String, Map<String, Double>> theoreticalValues = new HashMap<>(); // sample type (C/P/U) -> taxa -> value
		
		try (BufferedReader reader = new BufferedReader(new FileReader(paperDataFile))) {
			String line;
			boolean inRelativeTable = false;
			
			while ((line = reader.readLine()) != null) {
				line = line.trim();
				
				if (line.contains("All data Normalized to %")) {
					inRelativeTable = true;
					continue;
				}
				
				if (line.contains("Label\tSpecies")) {
					continue; // Skip header
				}
				
				if (line.isEmpty() || line.contains("Supplementary Data") || line.contains("basis for Figure") ||
					line.contains("Input Protein Amount")) {
					continue;
				}
				
				String[] parts = line.split("\t");
				if (parts.length < 8) {
					continue;
				}
				
				String label = parts[0].trim();
				if (label.isEmpty() || label.equals("Label")) {
					continue;
				}
				
				// Read relative values (normalized to %)
				// Column indices: Label=0, Species=1, Protein Amount: U=2, C=3, P=4
				if (inRelativeTable && parts.length > 4) {
					try {
						double u = parseDouble(parts[2].trim());
						double c = parseDouble(parts[3].trim());
						double p = parseDouble(parts[4].trim());
						
						if (!Double.isNaN(u) && u > 0) {
							theoreticalValues.computeIfAbsent("U", k -> new HashMap<>()).put(label, u);
						}
						if (!Double.isNaN(c) && c > 0) {
							theoreticalValues.computeIfAbsent("C", k -> new HashMap<>()).put(label, c);
						}
						if (!Double.isNaN(p) && p > 0) {
							theoreticalValues.computeIfAbsent("P", k -> new HashMap<>()).put(label, p);
						}
					} catch (NumberFormatException e) {
						// Skip
					}
				}
			}
		}
		
		System.out.println("Found theoretical values for " + theoreticalValues.size() + " sample types");
		
		// Step 3b: Filter and recalculate theoretical relative abundances using only the 19 specified taxa
		Set<String> validTaxa = new HashSet<>(Arrays.asList(
			"259", "841", "AK199", "CV", "VF", "BS", "KF7", "HB2", "PD", 
			"Am2", "NV", "BXL", "ATN", "Pfl", "CRH", "SMS", "K12", "LT2", "Cup"
		));
		
		System.out.println("Filtering theoretical values to " + validTaxa.size() + " specified taxa...");
		Map<String, Map<String, Double>> filteredTheoreticalValues = new HashMap<>(); // sample type -> taxa -> value
		
		// Filter and recalculate relative abundances for each sample type
		for (Map.Entry<String, Map<String, Double>> sampleTypeEntry : theoreticalValues.entrySet()) {
			String sampleType = sampleTypeEntry.getKey();
			Map<String, Double> originalValues = sampleTypeEntry.getValue();
			
			// Filter to only valid taxa
			Map<String, Double> filteredValues = new HashMap<>();
			double total = 0.0;
			for (Map.Entry<String, Double> entry : originalValues.entrySet()) {
				String taxa = entry.getKey();
				if (validTaxa.contains(taxa)) {
					double value = entry.getValue();
					filteredValues.put(taxa, value);
					total += value;
				}
			}
			
			// Recalculate relative abundances (normalize to 100%)
			if (total > 0) {
				Map<String, Double> normalizedValues = new HashMap<>();
				for (Map.Entry<String, Double> entry : filteredValues.entrySet()) {
					double normalized = (entry.getValue() / total) * 100.0;
					normalizedValues.put(entry.getKey(), normalized);
				}
				filteredTheoreticalValues.put(sampleType, normalizedValues);
				System.out.println("  Sample type " + sampleType + ": " + normalizedValues.size() + " taxa, total normalized to 100%");
			}
		}
		
		// Use filtered theoretical values for correlation and abundance table
		theoreticalValues = filteredTheoreticalValues;
		
		// Step 4: Write abundance tables for each method
		System.out.println("\nWriting abundance tables...");
		
		// Create output file for abundance tables (same directory, different name)
		File abundanceOutputFile = new File(outputFile.getParent(), 
			outputFile.getName().replace(".tsv", "_abundances.tsv").replace(".txt", "_abundances.tsv"));
		
		PrintWriter abundanceWriter = new PrintWriter(abundanceOutputFile);
		abundanceWriter.println("Method\tSample\tTaxa\tAbundance\tTheoretical\tPeptide_Count");
		
		// Collect all unique taxa across all methods and samples
		Set<String> allTaxa = new HashSet<>();
		for (String method : Arrays.asList("Sum", "Median", "Mean", "Max", "Sum_Distinct")) {
			Map<String, Map<String, Double>> methodAbundances = abundances.get(method);
			for (Map<String, Double> sampleAbundances : methodAbundances.values()) {
				allTaxa.addAll(sampleAbundances.keySet());
			}
		}
		
		// Write abundance data for each method
		for (String method : Arrays.asList("Sum", "Median", "Mean", "Max", "Sum_Distinct")) {
			Map<String, Map<String, Double>> methodAbundances = abundances.get(method);
			Map<String, Map<String, Integer>> methodPeptideCounts = peptideCounts.get(method);
			
			// Sort samples for consistent output
			ArrayList<String> sortedSamples = new ArrayList<>(methodAbundances.keySet());
			Collections.sort(sortedSamples);
			
			for (String sample : sortedSamples) {
				Map<String, Double> sampleAbundances = methodAbundances.get(sample);
				Map<String, Integer> samplePeptideCounts = methodPeptideCounts.get(sample);
				
				// Get sample type for theoretical value lookup
				String sampleType = extractSampleType(sample); // C, P, or U
				Map<String, Double> theoreticalForSampleType = null;
				if (sampleType != null && theoreticalValues.containsKey(sampleType)) {
					theoreticalForSampleType = theoreticalValues.get(sampleType);
				}
				
				// Sort taxa for consistent output
				ArrayList<String> sortedTaxa = new ArrayList<>(sampleAbundances.keySet());
				Collections.sort(sortedTaxa);
				
				for (String taxa : sortedTaxa) {
					double abundance = sampleAbundances.get(taxa);
					int peptideCount = (samplePeptideCounts != null && samplePeptideCounts.containsKey(taxa)) 
						? samplePeptideCounts.get(taxa) : 0;
					
					// Look up theoretical value
					double theoretical = Double.NaN;
					if (theoreticalForSampleType != null && theoreticalForSampleType.containsKey(taxa)) {
						theoretical = theoreticalForSampleType.get(taxa);
					}
					
					if (!Double.isNaN(theoretical)) {
						abundanceWriter.printf("%s\t%s\t%s\t%.6f\t%.6f\t%d\n", 
							method, sample, taxa, abundance, theoretical, peptideCount);
					} else {
						abundanceWriter.printf("%s\t%s\t%s\t%.6f\tNA\t%d\n", 
							method, sample, taxa, abundance, peptideCount);
					}
				}
			}
		}
		
		abundanceWriter.close();
		System.out.println("Abundance tables written to: " + abundanceOutputFile.getPath());
		
		// Step 5: Calculate correlations
		System.out.println("\nCalculating correlations...");
		
		PrintWriter writer = new PrintWriter(outputFile);
		writer.println("Method\tSample_Type\tSample\tCorrelation\tN");
		
		// Map samples to theoretical types: C1-C4 -> C, P1-P4 -> P, U1-U4 -> U
		for (String method : Arrays.asList("Sum", "Median", "Mean", "Max", "Sum_Distinct")) {
			Map<String, Map<String, Double>> methodAbundances = abundances.get(method);
			
			for (String sample : methodAbundances.keySet()) {
				String sampleType = extractSampleType(sample); // C, P, or U
				if (sampleType == null || !theoreticalValues.containsKey(sampleType)) {
					continue;
				}
				
				Map<String, Double> experimental = methodAbundances.get(sample);
				Map<String, Double> theoretical = theoreticalValues.get(sampleType);
				
				// Collect matching taxa
				ArrayList<Double> expValues = new ArrayList<>();
				ArrayList<Double> theoValues = new ArrayList<>();
				
				for (String taxa : experimental.keySet()) {
					if (theoretical.containsKey(taxa)) {
						expValues.add(experimental.get(taxa));
						theoValues.add(theoretical.get(taxa));
					}
				}
				
				if (expValues.size() >= 2) {
					double correlation = calculatePearsonCorrelation(theoValues, expValues);
					if (!Double.isNaN(correlation)) {
						writer.printf("%s\t%s\t%s\t%.6f\t%d\n", method, sampleType, sample, correlation, expValues.size());
					}
				}
			}
		}
		
		writer.close();
		System.out.println("Correlation results written to: " + outputFile.getPath());
	}
	
	/**
	 * Extract sample name from column header (e.g., "Run1_C1_2000ng.mzML" -> "C1")
	 */
	private String extractSampleName(String header) {
		// Pattern: Run1_C1_2000ng.mzML -> C1
		// Pattern: Run1_P1_96ng.mzML -> P1
		// Pattern: Run1_U1_2000ng.mzML -> U1
		if (header.contains("Run1_")) {
			int start = header.indexOf("Run1_") + 5;
			int end = header.indexOf("_", start);
			if (end < 0) {
				end = header.indexOf(".", start);
			}
			if (end > start) {
				return header.substring(start, end);
			}
		}
		return null;
	}
	
	/**
	 * Extract sample type from sample name (e.g., "C1" -> "C", "P2" -> "P", "U3" -> "U")
	 */
	private String extractSampleType(String sample) {
		if (sample != null && sample.length() > 0) {
			char firstChar = sample.charAt(0);
			if (firstChar == 'C' || firstChar == 'P' || firstChar == 'U') {
				return String.valueOf(firstChar);
			}
		}
		return null;
	}
	
	/**
	 * Calculate correlation between theoretical and experimental values from paper_data.txt
	 * Theoretical values: Input Protein Amount (U, C, P) and Input Cell Number (U, C, P)
	 * Experimental values: PSMs, Intensities (Unique Peptides), Intensities (Razor+Unique Peptides)
	 * Mapping: U (theoretical) matches U1, U2, U3, U4 (experimental), etc.
	 * Correlation is calculated across all species for each theoretical-experimental comparison
	 */
	public void calculateCorrelationFromPaperData(File paperDataFile, File outputFile) throws IOException {
		if (!paperDataFile.exists() || !paperDataFile.isFile()) {
			throw new IOException("Paper data file does not exist: " + paperDataFile.getAbsolutePath());
		}
		
		System.out.println("=== Correlation Analysis from Paper Data ===");
		System.out.println("Reading file: " + paperDataFile.getPath());
		
		// Structure to store all data for correlation calculation
		// Map: tableType -> theoreticalType -> method -> theoretical values list, experimental values list
		Map<String, Map<String, Map<String, CorrelationData>>> correlationData = new HashMap<>();
		
		// Read the file
		try (BufferedReader reader = new BufferedReader(new FileReader(paperDataFile))) {
			String line;
			String tableType = null; // "Absolute" or "Relative"
			
			while ((line = reader.readLine()) != null) {
				line = line.trim();
				
				if (line.isEmpty()) {
					continue;
				}
				
				// Detect table type
				if (line.contains("All data Normalized to %") || line.contains("Normalized to %")) {
					tableType = "Relative";
					continue;
				}
				
				// Skip header lines
				if (line.contains("Supplementary Data") || line.contains("basis for Figure") || 
					line.contains("Input Protein Amount") || line.contains("Label\tSpecies")) {
					if (tableType == null && line.contains("Input Protein Amount")) {
						tableType = "Absolute";
					}
					continue;
				}
				
				// Parse data rows
				String[] parts = line.split("\t");
				if (parts.length < 44) { // Need at least 44 columns
					continue;
				}
				
				String label = parts[0].trim();
				if (label.isEmpty() || label.equals("Label")) {
					continue; // Skip header row
				}
				
				// Extract theoretical and experimental values for correlation
				// Column indices: Label=0, Species=1, Protein Amount: U=2, C=3, P=4, Cell Number: U=5, C=6, P=7
				// Experimental: PSMs C1=8-11, P1=12-15, U1=16-19
				// Unique: C1=20-23, P1=24-27, U1=28-31
				// Razor+Unique: C1=32-35, P1=36-39, U1=40-43
				
				if (tableType != null) {
					// Collect data for PSMs (start at column 8)
					collectCorrelationData(correlationData, tableType, parts, 8, "PSMs");
					
					// Collect data for Intensities (Unique Peptides) (start at column 20)
					collectCorrelationData(correlationData, tableType, parts, 20, "Intensities_Unique");
					
					// Collect data for Intensities (Razor+Unique Peptides) (start at column 32)
					collectCorrelationData(correlationData, tableType, parts, 32, "Intensities_RazorUnique");
				}
			}
		}
		
		// Calculate and write correlations
		PrintWriter writer = new PrintWriter(outputFile);
		writer.println("Table_Type\tTheoretical_Type\tExperimental_Type\tCorrelation\tN");
		
		for (Map.Entry<String, Map<String, Map<String, CorrelationData>>> tableEntry : correlationData.entrySet()) {
			String tableType = tableEntry.getKey();
			Map<String, Map<String, CorrelationData>> theoreticalMap = tableEntry.getValue();
			
			for (Map.Entry<String, Map<String, CorrelationData>> theoreticalEntry : theoreticalMap.entrySet()) {
				String theoreticalType = theoreticalEntry.getKey();
				Map<String, CorrelationData> methodMap = theoreticalEntry.getValue();
				
				for (Map.Entry<String, CorrelationData> methodEntry : methodMap.entrySet()) {
					String methodName = methodEntry.getKey();
					CorrelationData data = methodEntry.getValue();
					
					if (data.theoretical.size() >= 2 && data.experimental.size() >= 2 && 
						data.theoretical.size() == data.experimental.size()) {
						double correlation = calculatePearsonCorrelation(data.theoretical, data.experimental);
						if (!Double.isNaN(correlation)) {
							writer.printf("%s\t%s\t%s\t%.6f\t%d\n", 
								tableType, theoreticalType, methodName, correlation, data.theoretical.size());
						}
					}
				}
			}
		}
		
		writer.close();
		System.out.println("Correlation results written to: " + outputFile.getPath());
	}
	
	/**
	 * Helper class to store theoretical and experimental values for correlation
	 */
	private static class CorrelationData {
		ArrayList<Double> theoretical = new ArrayList<>();
		ArrayList<Double> experimental = new ArrayList<>();
	}
	
	/**
	 * Collect correlation data for a specific experimental method
	 */
	private void collectCorrelationData(Map<String, Map<String, Map<String, CorrelationData>>> correlationData,
			String tableType, String[] parts, int startCol, String methodName) {
		
		// Initialize structure
		correlationData.computeIfAbsent(tableType, k -> new HashMap<>());
		Map<String, Map<String, CorrelationData>> theoreticalMap = correlationData.get(tableType);
		
		// Column indices: Protein Amount: U=2, C=3, P=4, Cell Number: U=5, C=6, P=7
		// Experimental: C1=startCol, C2=startCol+1, C3=startCol+2, C4=startCol+3
		//                P1=startCol+4, P2=startCol+5, P3=startCol+6, P4=startCol+7
		//                U1=startCol+8, U2=startCol+9, U3=startCol+10, U4=startCol+11
		
		// For U (uneven): theoretical column 2 (Protein) or 5 (Cell), experimental U1-U4 (startCol+8 to startCol+11)
		if (parts.length > startCol + 11) {
			// Protein Amount U
			double theoreticalU_Protein = parseDouble(parts[2].trim());
			if (!Double.isNaN(theoreticalU_Protein) && theoreticalU_Protein > 0) {
				ArrayList<Double> experimental = extractExperimentalValues(parts, startCol + 8, 4);
				if (!experimental.isEmpty()) {
					theoreticalMap.computeIfAbsent("Protein_Amount_U", k -> new HashMap<>())
						.computeIfAbsent(methodName, k -> new CorrelationData());
					CorrelationData data = theoreticalMap.get("Protein_Amount_U").get(methodName);
					for (double exp : experimental) {
						data.theoretical.add(theoreticalU_Protein);
						data.experimental.add(exp);
					}
				}
			}
			
			// Cell Number U
			double theoreticalU_Cell = parseDouble(parts[5].trim());
			if (!Double.isNaN(theoreticalU_Cell) && theoreticalU_Cell > 0) {
				ArrayList<Double> experimental = extractExperimentalValues(parts, startCol + 8, 4);
				if (!experimental.isEmpty()) {
					theoreticalMap.computeIfAbsent("Cell_Number_U", k -> new HashMap<>())
						.computeIfAbsent(methodName, k -> new CorrelationData());
					CorrelationData data = theoreticalMap.get("Cell_Number_U").get(methodName);
					for (double exp : experimental) {
						data.theoretical.add(theoreticalU_Cell);
						data.experimental.add(exp);
					}
				}
			}
		}
		
		// For C (equal cell number): theoretical column 3 (Protein) or 6 (Cell), experimental C1-C4 (startCol to startCol+3)
		if (parts.length > startCol + 3) {
			// Protein Amount C
			double theoreticalC_Protein = parseDouble(parts[3].trim());
			if (!Double.isNaN(theoreticalC_Protein) && theoreticalC_Protein > 0) {
				ArrayList<Double> experimental = extractExperimentalValues(parts, startCol, 4);
				if (!experimental.isEmpty()) {
					theoreticalMap.computeIfAbsent("Protein_Amount_C", k -> new HashMap<>())
						.computeIfAbsent(methodName, k -> new CorrelationData());
					CorrelationData data = theoreticalMap.get("Protein_Amount_C").get(methodName);
					for (double exp : experimental) {
						data.theoretical.add(theoreticalC_Protein);
						data.experimental.add(exp);
					}
				}
			}
			
			// Cell Number C
			double theoreticalC_Cell = parseDouble(parts[6].trim());
			if (!Double.isNaN(theoreticalC_Cell) && theoreticalC_Cell > 0) {
				ArrayList<Double> experimental = extractExperimentalValues(parts, startCol, 4);
				if (!experimental.isEmpty()) {
					theoreticalMap.computeIfAbsent("Cell_Number_C", k -> new HashMap<>())
						.computeIfAbsent(methodName, k -> new CorrelationData());
					CorrelationData data = theoreticalMap.get("Cell_Number_C").get(methodName);
					for (double exp : experimental) {
						data.theoretical.add(theoreticalC_Cell);
						data.experimental.add(exp);
					}
				}
			}
		}
		
		// For P (equal protein amounts): theoretical column 4 (Protein) or 7 (Cell), experimental P1-P4 (startCol+4 to startCol+7)
		if (parts.length > startCol + 7) {
			// Protein Amount P
			double theoreticalP_Protein = parseDouble(parts[4].trim());
			if (!Double.isNaN(theoreticalP_Protein) && theoreticalP_Protein > 0) {
				ArrayList<Double> experimental = extractExperimentalValues(parts, startCol + 4, 4);
				if (!experimental.isEmpty()) {
					theoreticalMap.computeIfAbsent("Protein_Amount_P", k -> new HashMap<>())
						.computeIfAbsent(methodName, k -> new CorrelationData());
					CorrelationData data = theoreticalMap.get("Protein_Amount_P").get(methodName);
					for (double exp : experimental) {
						data.theoretical.add(theoreticalP_Protein);
						data.experimental.add(exp);
					}
				}
			}
			
			// Cell Number P
			double theoreticalP_Cell = parseDouble(parts[7].trim());
			if (!Double.isNaN(theoreticalP_Cell) && theoreticalP_Cell > 0) {
				ArrayList<Double> experimental = extractExperimentalValues(parts, startCol + 4, 4);
				if (!experimental.isEmpty()) {
					theoreticalMap.computeIfAbsent("Cell_Number_P", k -> new HashMap<>())
						.computeIfAbsent(methodName, k -> new CorrelationData());
					CorrelationData data = theoreticalMap.get("Cell_Number_P").get(methodName);
					for (double exp : experimental) {
						data.theoretical.add(theoreticalP_Cell);
						data.experimental.add(exp);
					}
				}
			}
		}
	}
	
	/**
	 * Extract experimental values from specified columns
	 */
	private ArrayList<Double> extractExperimentalValues(String[] parts, int startCol, int count) {
		ArrayList<Double> values = new ArrayList<>();
		for (int i = startCol; i < startCol + count && i < parts.length; i++) {
			double val = parseDouble(parts[i].trim());
			if (!Double.isNaN(val) && val >= 0) {
				values.add(val);
			}
		}
		return values;
	}
	
	/**
	 * Parse double value, handling scientific notation (e.g., "5.79E+09")
	 */
	private double parseDouble(String str) {
		if (str == null || str.trim().isEmpty() || str.equals("NA") || str.equals("N/A")) {
			return Double.NaN;
		}
		try {
			return Double.parseDouble(str.trim());
		} catch (NumberFormatException e) {
			return Double.NaN;
		}
	}
	
	/**
	 * Calculate Pearson correlation coefficient between two lists
	 */
	private double calculatePearsonCorrelation(ArrayList<Double> x, ArrayList<Double> y) {
		if (x.size() != y.size() || x.size() < 2) {
			return Double.NaN;
		}
		
		// Calculate means
		double meanX = 0.0;
		double meanY = 0.0;
		for (int i = 0; i < x.size(); i++) {
			meanX += x.get(i);
			meanY += y.get(i);
		}
		meanX /= x.size();
		meanY /= y.size();
		
		// Calculate numerator and denominators
		double numerator = 0.0;
		double sumSqX = 0.0;
		double sumSqY = 0.0;
		
		for (int i = 0; i < x.size(); i++) {
			double diffX = x.get(i) - meanX;
			double diffY = y.get(i) - meanY;
			numerator += diffX * diffY;
			sumSqX += diffX * diffX;
			sumSqY += diffY * diffY;
		}
		
		double denominator = Math.sqrt(sumSqX * sumSqY);
		if (denominator == 0.0) {
			return Double.NaN;
		}
		
		return numerator / denominator;
	}
	
	/**
	 * Write peptide detectability statistics
	 * Outputs a table with: Peptide, Detectability (from predicted files), Category (marker type or Others)
	 */
	private void writePeptideDetectability(Map<String, PeptideData> peptideData, 
			Map<String, Double> peptideToDetectability, File faaFolder, File outputFile) throws IOException {
		
		// Load protein-to-category mapping
		Map<String, String> proteinToCategory = loadMarkerProteinCategories(faaFolder);
		
		PrintWriter writer = new PrintWriter(outputFile);
		writer.println("Peptide\tDetectability\tCategory");
		
		// Process each peptide
		for (Map.Entry<String, PeptideData> entry : peptideData.entrySet()) {
			String peptide = entry.getKey();
			PeptideData data = entry.getValue();
			
			// Get detectability (from predicted files)
			Double detectability = peptideToDetectability.get(peptide);
			if (detectability != null) {
				// Determine category based on proteins
				String category = determinePeptideCategory(data.proteinIds, proteinToCategory);
				writer.printf("%s\t%.6f\t%s\n", peptide, detectability, category);
			}			
		}
		
		writer.close();
	}
	
	/**
	 * Determine the category for a peptide based on its proteins
	 * Returns: "Anchor", "HighRes", "Growth", "Housekeeping", or "Others"
	 */
	private String determinePeptideCategory(Set<String> proteinIds, Map<String, String> proteinToCategory) {
		// Count categories for this peptide's proteins
		Map<String, Integer> categoryCounts = new HashMap<>();
		categoryCounts.put("Anchor", 0);
		categoryCounts.put("HighRes", 0);
		categoryCounts.put("Growth", 0);
		categoryCounts.put("Housekeeping", 0);
		categoryCounts.put("Others", 0);
		
		for (String proteinId : proteinIds) {
			String category = proteinToCategory.get(proteinId);
			if (category != null) {
				categoryCounts.put(category, categoryCounts.getOrDefault(category, 0) + 1);
			} else {
				categoryCounts.put("Others", categoryCounts.get("Others") + 1);
			}
		}
		
		// Find the most common category (excluding Others if there are any marker categories)
		int maxCount = 0;
		String mostCommonCategory = "Others";
		
		// First, check marker categories
		for (String markerCategory : Arrays.asList("Anchor", "HighRes", "Growth", "Housekeeping")) {
			int count = categoryCounts.get(markerCategory);
			if (count > maxCount) {
				maxCount = count;
				mostCommonCategory = markerCategory;
			}
		}
		
		// If no marker proteins found, use Others
		if (maxCount == 0) {
			mostCommonCategory = "Others";
		}
		
		return mostCommonCategory;
	}
	
	public void getGenomeCount(String in) {
		File[] files = (new File(in)).listFiles();
		for (int i = 0; i < files.length; i++) {
			try (BufferedReader reader = new BufferedReader(new FileReader(files[i]))) {
				String line = reader.readLine();
				HashSet<String> set1 = new HashSet<String>();
				HashSet<String> set2 = new HashSet<String>();
				while((line=reader.readLine())!=null) {
					String[] cs = line.split("\t");
					set1.add(cs[0]);
					set2.add(cs[1]);
				}
				reader.close();
				System.out.println(files[i].getName()+"\t"+set1.size()+"\t"+set2.size());
			} catch (IOException e) {

			}
		}
	}
	
	/**
	 * Test method to demonstrate protein information extraction
	 */
	public static void testExtractProteinInfo() {
		DatasetMarkerAnalysis analyzer = new DatasetMarkerAnalysis();
		
		// Test case 1: OS= format
		String header1 = ">ATN_Q8U9B0 RGMG3_AGRFC Putative ribose/galactose/methyl galactoside import ATP-binding protein 3 OS=Agrobacterium fabrum (strain C58 / ATCC 33970) GN=Atu3818 PE=3 SV=2";
		ProteinInfo info1 = analyzer.extractProteinInfo(header1);
		System.out.println("Test 1:");
		System.out.println("  Input: " + header1);
		System.out.println("  Output: " + info1);
		System.out.println();
		
		// Test case 2: Bracket format
		String header2 = ">CRH_XP_001690518.1 hypothetical protein CHLREDRAFT_169266 [Chlamydomonas reinhardtii]";
		ProteinInfo info2 = analyzer.extractProteinInfo(header2);
		System.out.println("Test 2:");
		System.out.println("  Input: " + header2);
		System.out.println("  Output: " + info2);
		System.out.println();
	}
	
	public void exportPeptideFasta(String in, String out, double threshold) {

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
	
	public void getSimpleTaxaInfo(String in, String out) {
		try (BufferedReader reader = new BufferedReader(new FileReader(in))) {
			PrintWriter writer = new PrintWriter(out);
			String line = reader.readLine();
			HashMap<String, String> map  = new HashMap<String, String>();
			while((line=reader.readLine())!=null) {
				String[] cs = line.split("\t");
				if(!map.containsKey(cs[13])) {
					map.put(cs[13], cs[14]);
				}
			}
			reader.close();
			for(String key: map.keySet()) {
				writer.println(key+"\t"+map.get(key));
			}
			writer.close();
		} catch (IOException e) {

		}
	}

	public void evaluateGenomeSelection(File genomesFile, File hapGenomesFile) {
		Map<String, Set<String>> taxonToSpeciesIds = new HashMap<>();
		Map<String, Set<String>> taxonToGenusIds = new HashMap<>();
		Set<String> allSpeciesIds = new HashSet<>();
		Set<String> allGenusIds = new HashSet<>();

		try (BufferedReader reader = new BufferedReader(new FileReader(genomesFile))) {
			String line = reader.readLine(); // header
			while ((line = reader.readLine()) != null) {
				if (line.trim().isEmpty()) {
					continue;
				}
				String[] parts = line.split(",", -1);
				if (parts.length < 4) {
					continue;
				}
				String taxon = parts[1].trim();
				Set<String> speciesIds = parseGenomeIdList(parts[2]);
				Set<String> genusIds = parseGenomeIdList(parts[3]);

				taxonToSpeciesIds.put(taxon, speciesIds);
				taxonToGenusIds.put(taxon, genusIds);
				allSpeciesIds.addAll(speciesIds);
				allGenusIds.addAll(genusIds);
			}
		} catch (IOException e) {
			System.out.println("Error reading genomes.txt: " + e.getMessage());
			return;
		}

		class SelectedGenome {
			final int rank;
			final String genomeId;
			final int increment;

			SelectedGenome(int rank, String genomeId, int increment) {
				this.rank = rank;
				this.genomeId = genomeId;
				this.increment = increment;
			}
		}

		Set<String> selectedGenomes = new HashSet<>();
		ArrayList<SelectedGenome> selectedGenomeList = new ArrayList<>();
		try (BufferedReader reader = new BufferedReader(new FileReader(hapGenomesFile))) {
			String line = reader.readLine(); // header
			while ((line = reader.readLine()) != null) {
				if (line.trim().isEmpty()) {
					continue;
				}
				String[] parts = line.split("\t");
				if (parts.length < 6) {
					continue;
				}
				int rank;
				try {
					rank = Integer.parseInt(parts[0].trim());
				} catch (NumberFormatException e) {
					continue;
				}
				String genomeId = parts[1].trim();
				String incrementStr = parts[5].trim();
				int increment = 0;
				try {
					increment = Integer.parseInt(incrementStr);
				} catch (NumberFormatException e) {
					// Leave as 0
				}
				selectedGenomes.add(genomeId);
				selectedGenomeList.add(new SelectedGenome(rank, genomeId, increment));
				if (increment == 1) {
					break;
				}
			}
		} catch (IOException e) {
			System.out.println("Error reading hap_genomes.tsv: " + e.getMessage());
			return;
		}

		Set<String> speciesMatchedGenomes = new HashSet<>();
		Set<String> genusMatchedGenomes = new HashSet<>();
		for (String genomeId : selectedGenomes) {
			if (allSpeciesIds.contains(genomeId)) {
				speciesMatchedGenomes.add(genomeId);
			}
			if (allGenusIds.contains(genomeId)) {
				genusMatchedGenomes.add(genomeId);
			}
		}

		int speciesMatchCount = speciesMatchedGenomes.size();
		int genusMatchCount = genusMatchedGenomes.size();
		Set<String> bothMatchedGenomes = new HashSet<>(speciesMatchedGenomes);
		bothMatchedGenomes.retainAll(genusMatchedGenomes);
		Set<String> genusOnlyGenomes = new HashSet<>(genusMatchedGenomes);
		genusOnlyGenomes.removeAll(speciesMatchedGenomes);
		Set<String> incorrectMatchedGenomes = new HashSet<>(selectedGenomes);
		incorrectMatchedGenomes.removeAll(speciesMatchedGenomes);
		incorrectMatchedGenomes.removeAll(genusMatchedGenomes);

		int notFoundTaxaCount = 0;
		ArrayList<String> notFoundTaxa = new ArrayList<>();
		for (String taxon : taxonToSpeciesIds.keySet()) {
			Set<String> speciesIds = taxonToSpeciesIds.getOrDefault(taxon, Collections.emptySet());
			Set<String> genusIds = taxonToGenusIds.getOrDefault(taxon, Collections.emptySet());
			boolean found = false;
			for (String genomeId : selectedGenomes) {
				if (speciesIds.contains(genomeId) || genusIds.contains(genomeId)) {
					found = true;
					break;
				}
			}
			if (!found) {
				notFoundTaxaCount++;
				notFoundTaxa.add(taxon);
			}
		}

		Collections.sort(notFoundTaxa);
		System.out.println("Selected genomes (up to first Increment==1): " + selectedGenomes.size());
		System.out.println("Species match genomes: " + speciesMatchCount);
		System.out.println("Genus match genomes: " + genusMatchCount);
		System.out.println("Genus-only genomes: " + genusOnlyGenomes.size());
		System.out.println("Species+Genus genomes: " + bothMatchedGenomes.size());
		if (!genusMatchedGenomes.isEmpty()) {
			ArrayList<String> sortedGenusMatched = new ArrayList<>(genusMatchedGenomes);
			Collections.sort(sortedGenusMatched);
			System.out.println("Genus match genome IDs: " + String.join(", ", sortedGenusMatched));
		}
		if (!speciesMatchedGenomes.isEmpty()) {
			ArrayList<String> sortedSpeciesMatched = new ArrayList<>(speciesMatchedGenomes);
			Collections.sort(sortedSpeciesMatched);
			System.out.println("Species match genome IDs: " + String.join(", ", sortedSpeciesMatched));
		}
		System.out.println("Incorrect matched genomes: " + incorrectMatchedGenomes.size());
		if (!incorrectMatchedGenomes.isEmpty()) {
			ArrayList<String> sortedIncorrect = new ArrayList<>(incorrectMatchedGenomes);
			Collections.sort(sortedIncorrect);
			System.out.println("Incorrect genome IDs: " + String.join(", ", sortedIncorrect));
		}
		System.out.println("Taxa not found: " + notFoundTaxaCount);
		if (!notFoundTaxa.isEmpty()) {
			System.out.println("Not found taxa: " + String.join(", ", notFoundTaxa));
		}
		if (!selectedGenomeList.isEmpty()) {
			Collections.sort(selectedGenomeList, Comparator.comparingInt(a -> a.rank));
			System.out.println("Rank\tGenomeID\tStatus");
			for (SelectedGenome entry : selectedGenomeList) {
				String status;
				if (speciesMatchedGenomes.contains(entry.genomeId)) {
					status = "Species";
				} else if (genusMatchedGenomes.contains(entry.genomeId)) {
					status = "Genus";
				} else {
					status = "Incorrect";
				}
				System.out.println(entry.rank + "\t" + entry.genomeId + "\t" + status);
			}
		}
	}

	private Set<String> parseGenomeIdList(String raw) {
		Set<String> ids = new HashSet<>();
		String value = raw == null ? "" : raw.trim();
		if (value.isEmpty() || "None".equalsIgnoreCase(value)) {
			return ids;
		}
		String[] parts = value.split(";");
		for (String part : parts) {
			String id = part.trim();
			if (!id.isEmpty() && !"None".equalsIgnoreCase(id)) {
				ids.add(id);
			}
		}
		return ids;
	}
	
	public void peptideTest(String in, File markerFasta, String output) {
		int[] counts = new int[10];
		HashSet<String> pepset = new HashSet<String>();
		Map<String, String> proteinToCategory = new HashMap<>();
		try {
			proteinToCategory = loadMarkerProteinCategoriesFromFasta(markerFasta);
		} catch (IOException e) {
			// Ignore and continue with empty mapping
		}
		System.out.println("proteinToCategory\t" + proteinToCategory.size());
		Map<String, Map<String, Set<String>>> genomeCategoryPeptides = new HashMap<>();
		Map<String, Map<String, Set<String>>> genomeCategoryUniquePeptides = new HashMap<>();
		Map<String, Set<String>> peptideGenomes = new HashMap<>();
		Map<String, Map<String, Set<String>>> peptideGenomeCategories = new HashMap<>();
		try (BufferedReader hapReader = new BufferedReader(new FileReader(in))) {
			String headerLine = hapReader.readLine(); // Read header
			String[] header = headerLine.split("\t");

			// Find the index of "proteins" column
			int peptidesIndex = -1;
			int proteinsIndex = -1;
			int labelIndex = -1;
			int peptideQIndex = -1;
			for (int j = 0; j < header.length; j++) {
				if (header[j].toLowerCase().equals("peptide")) {
					peptidesIndex = j;
				} else if (header[j].toLowerCase().equals("proteins")) {
					proteinsIndex = j;
				} else if (header[j].toLowerCase().equals("label")) {
					labelIndex = j;
				} else if (header[j].toLowerCase().equals("peptide_q")) {
					peptideQIndex = j;
				}
			}

			String line;
			while ((line = hapReader.readLine()) != null) {
				String[] cs = line.split("\t");
				if (cs.length > labelIndex) {
					if (!cs[labelIndex].equals("1")) {
						continue;
					}
				} else {
					continue;
				}

				if (cs.length > peptideQIndex) {
					if (Double.parseDouble(cs[peptideQIndex]) > 0.01) {
						continue;
					}
				} else {
					continue;
				}

				if (cs.length > proteinsIndex) {
					String peptide = cs[peptidesIndex];
					if (!pepset.contains(peptide)) {
						String[] pros = cs[proteinsIndex].split(";");
						if (pros.length >= counts.length) {
							counts[counts.length - 1]++;
						} else {
							counts[pros.length]++;
						}
						pepset.add(peptide);

						for (String pro : pros) {
							String proteinId = pro.trim();
							if (proteinId.isEmpty()) {
								continue;
							}
							proteinId = proteinId.substring(0, proteinId.lastIndexOf("_"));
							String category = proteinToCategory.get(proteinId);
							if (category == null) {
								continue;
							}
							String genomeId = extractGenomeId(proteinId);

							peptideGenomes.computeIfAbsent(peptide, k -> new HashSet<>()).add(genomeId);
							Map<String, Set<String>> genomeMap = genomeCategoryPeptides.computeIfAbsent(genomeId,
									k -> new HashMap<>());
							genomeMap.computeIfAbsent(category, k -> new HashSet<>()).add(peptide);

							Map<String, Set<String>> genomeCategories = peptideGenomeCategories.computeIfAbsent(peptide,
									k -> new HashMap<>());
							genomeCategories.computeIfAbsent(genomeId, k -> new HashSet<>()).add(category);
						}
					}
				}
			}

			System.out.println(pepset.size() + "\t" + Arrays.toString(counts));
			for (int i = 0; i < counts.length; i++) {
				System.out.print((double) counts[i] / (double) pepset.size() + "\t");
			}
			System.out.println();

			for (Map.Entry<String, Set<String>> entry : peptideGenomes.entrySet()) {
				Set<String> genomes = entry.getValue();
				if (genomes.size() != 1) {
					continue;
				}
				String peptide = entry.getKey();
				String genomeId = genomes.iterator().next();
				Map<String, Set<String>> genomeCategories = peptideGenomeCategories.get(peptide);
				if (genomeCategories == null) {
					continue;
				}
				Set<String> categories = genomeCategories.get(genomeId);
				if (categories == null) {
					continue;
				}
				Map<String, Set<String>> genomeUniqueMap = genomeCategoryUniquePeptides.computeIfAbsent(genomeId,
						k -> new HashMap<>());
				for (String category : categories) {
					genomeUniqueMap.computeIfAbsent(category, k -> new HashSet<>()).add(peptide);
				}
			}

			ArrayList<String> genomes = new ArrayList<>(genomeCategoryPeptides.keySet());
			Collections.sort(genomes);
			ArrayList<String> categories = new ArrayList<>(
					Arrays.asList("Anchor", "HighRes", "Growth", "Housekeeping"));

			PrintWriter writer = new PrintWriter(output);
			writer.println(
					"Genome\tAnchor\tHighRes\tGrowth\tHousekeeping\tUniqueAnchor\tUniqueHighRes\tUniqueGrowth\tUniqueHousekeeping\tTotalMarkers\tTotalUniqueMarkers");
			for (String genomeId : genomes) {
				Map<String, Set<String>> categoryPeptides = genomeCategoryPeptides.get(genomeId);
				Map<String, Set<String>> uniqueMap = genomeCategoryUniquePeptides.get(genomeId);
				int totalMarkers = 0;
				int totalUniqueMarkers = 0;
				StringBuilder row = new StringBuilder();
				row.append(genomeId);
				for (String category : categories) {
					int peptideCount = 0;
					if (categoryPeptides != null && categoryPeptides.containsKey(category)) {
						peptideCount = categoryPeptides.get(category).size();
					}
					totalMarkers += peptideCount;
					row.append("\t").append(peptideCount);
				}
				for (String category : categories) {
					int uniqueCount = 0;
					if (uniqueMap != null && uniqueMap.containsKey(category)) {
						uniqueCount = uniqueMap.get(category).size();
					}
					totalUniqueMarkers += uniqueCount;
					row.append("\t").append(uniqueCount);
				}
				row.append("\t").append(totalMarkers).append("\t").append(totalUniqueMarkers);
				writer.println(row.toString());
			}
			writer.close();

		} catch (IOException e) {

		}
	}
/*
	public void testMarkerGenomeSelectorEdgeCases(File markerFasta, File peptideTsv, File taxaTsv,
			File evidenceTsv, File rankedTsv, File selectedTsv, File diagnosticsTsv) throws IOException {
		runMarkerGenomeSelector(markerFasta, peptideTsv, taxaTsv, evidenceTsv, rankedTsv, selectedTsv);

		Map<String, String> proteinToCategory = loadMarkerProteinCategoriesFromFasta(markerFasta);
		Set<String> allMarkerPeptides = new HashSet<>();
		Map<String, Set<String>> peptideToGenomes = new HashMap<>();
		Map<String, Map<String, Set<String>>> genomeCategoryPeptides = new HashMap<>();
		Map<String, Map<String, Set<String>>> peptideUnionGenomeCategories = new HashMap<>();
		Map<String, Map<String, Set<String>>> peptideFirstGenomeCategories = new HashMap<>();
		Set<String> peptideFirstMarkerSeen = new HashSet<>();
		Set<String> noUnderscoreProteins = new HashSet<>();
		Set<String> noUnderscoreMarkerProteins = new HashSet<>();
		int totalLines = 0;
		int keptLines = 0;

		try (BufferedReader hapReader = new BufferedReader(new FileReader(peptideTsv))) {
			String headerLine = hapReader.readLine();
			if (headerLine == null) {
				throw new IOException("Empty peptide TSV: " + peptideTsv.getAbsolutePath());
			}
			String[] header = headerLine.split("\t");
			int peptidesIndex = -1;
			int proteinsIndex = -1;
			int labelIndex = -1;
			int peptideQIndex = -1;
			for (int j = 0; j < header.length; j++) {
				String h = header[j].toLowerCase();
				if (h.equals("peptide")) {
					peptidesIndex = j;
				} else if (h.equals("proteins")) {
					proteinsIndex = j;
				} else if (h.equals("label")) {
					labelIndex = j;
				} else if (h.equals("peptide_q")) {
					peptideQIndex = j;
				}
			}
			if (peptidesIndex < 0 || proteinsIndex < 0 || labelIndex < 0 || peptideQIndex < 0) {
				throw new IOException("Missing required columns: peptide/proteins/label/peptide_q");
			}
			String line;
			while ((line = hapReader.readLine()) != null) {
				totalLines++;
				String[] cs = line.split("\t");
				if (cs.length <= Math.max(Math.max(peptidesIndex, proteinsIndex), Math.max(labelIndex, peptideQIndex))) {
					continue;
				}
				if (!"1".equals(cs[labelIndex])) {
					continue;
				}
				if (Double.parseDouble(cs[peptideQIndex]) > 0.01) {
					continue;
				}
				String peptide = cs[peptidesIndex];
				String[] pros = cs[proteinsIndex].split(";");
				Map<String, Set<String>> lineGenomeCategories = new HashMap<>();
				boolean lineHasMarker = false;

				for (String pro : pros) {
					String proteinId = pro.trim();
					if (proteinId.isEmpty()) {
						continue;
					}
					int last = proteinId.lastIndexOf("_");
					if (last <= 0) {
						noUnderscoreProteins.add(proteinId);
						if (proteinToCategory.containsKey(proteinId)) {
							noUnderscoreMarkerProteins.add(proteinId);
						}
						continue;
					}
					String baseId = proteinId.substring(0, last);
					String category = proteinToCategory.get(baseId);
					if (category == null) {
						continue;
					}
					String genomeId = extractGenomeId(baseId);
					lineHasMarker = true;
					lineGenomeCategories.computeIfAbsent(genomeId, k -> new HashSet<>()).add(category);
				}

				if (lineHasMarker) {
					keptLines++;
					Map<String, Set<String>> unionMap = peptideUnionGenomeCategories.computeIfAbsent(peptide,
							k -> new HashMap<>());
					for (Map.Entry<String, Set<String>> entry : lineGenomeCategories.entrySet()) {
						unionMap.computeIfAbsent(entry.getKey(), k -> new HashSet<>()).addAll(entry.getValue());
					}
					if (!peptideFirstMarkerSeen.contains(peptide)) {
						peptideFirstMarkerSeen.add(peptide);
						peptideFirstGenomeCategories.put(peptide, copyGenomeCategories(lineGenomeCategories));
					}
				}

				if (allMarkerPeptides.contains(peptide)) {
					continue;
				}
				if (lineHasMarker) {
					allMarkerPeptides.add(peptide);
					for (Map.Entry<String, Set<String>> entry : lineGenomeCategories.entrySet()) {
						String genomeId = entry.getKey();
						peptideToGenomes.computeIfAbsent(peptide, k -> new HashSet<>()).add(genomeId);
						for (String category : entry.getValue()) {
							genomeCategoryPeptides.computeIfAbsent(genomeId, k -> new HashMap<>())
									.computeIfAbsent(category, k -> new HashSet<>())
									.add(peptide);
						}
					}
				}
			}
		}

		Set<String> hardFilteredGenomes = new HashSet<>();
		for (Map.Entry<String, Map<String, Set<String>>> entry : genomeCategoryPeptides.entrySet()) {
			Map<String, Set<String>> cat = entry.getValue();
			boolean passes = cat.containsKey("Anchor")
					&& (cat.containsKey("HighRes") || cat.containsKey("Growth"));
			if (passes) {
				hardFilteredGenomes.add(entry.getKey());
			}
		}
		Set<String> coverable = new HashSet<>();
		for (Map.Entry<String, Set<String>> entry : peptideToGenomes.entrySet()) {
			for (String genome : entry.getValue()) {
				if (hardFilteredGenomes.contains(genome)) {
					coverable.add(entry.getKey());
					break;
				}
			}
		}

		int lostPeptideCount = 0;
		ArrayList<String> lostExamples = new ArrayList<>();
		for (Map.Entry<String, Map<String, Set<String>>> entry : peptideUnionGenomeCategories.entrySet()) {
			String peptide = entry.getKey();
			Map<String, Set<String>> union = entry.getValue();
			Map<String, Set<String>> first = peptideFirstGenomeCategories.get(peptide);
			if (first == null) {
				continue;
			}
			boolean lost = false;
			for (Map.Entry<String, Set<String>> unionEntry : union.entrySet()) {
				Set<String> firstCats = first.get(unionEntry.getKey());
				if (firstCats == null || !firstCats.containsAll(unionEntry.getValue())) {
					lost = true;
					break;
				}
			}
			if (lost) {
				lostPeptideCount++;
				if (lostExamples.size() < 50) {
					lostExamples.add(peptide);
				}
			}
		}

		Map<String, Map<String, Set<String>>> genomeCategoryUniquePeptides = new HashMap<>();
		for (Map.Entry<String, Set<String>> entry : peptideToGenomes.entrySet()) {
			if (entry.getValue().size() != 1) {
				continue;
			}
			String peptide = entry.getKey();
			String genomeId = entry.getValue().iterator().next();
			Map<String, Set<String>> categories = peptideFirstGenomeCategories.get(peptide);
			if (categories == null) {
				continue;
			}
			Set<String> cats = categories.get(genomeId);
			if (cats == null) {
				continue;
			}
			Map<String, Set<String>> genomeUniqueMap = genomeCategoryUniquePeptides.computeIfAbsent(genomeId,
					k -> new HashMap<>());
			for (String category : cats) {
				genomeUniqueMap.computeIfAbsent(category, k -> new HashSet<>()).add(peptide);
			}
		}

		Set<String> filteredOutPeptides = new HashSet<>(allMarkerPeptides);
		filteredOutPeptides.removeAll(coverable);
		ArrayList<String> filteredOutList = new ArrayList<>(filteredOutPeptides);
		Collections.sort(filteredOutList);

		try (PrintWriter writer = new PrintWriter(diagnosticsTsv)) {
			writer.println("Metric\tValue");
			writer.println("Total_Lines\t" + totalLines);
			writer.println("Lines_With_Marker\t" + keptLines);
			writer.println("Total_Marker_Peptides\t" + allMarkerPeptides.size());
			writer.println("HardFiltered_Genomes\t" + hardFilteredGenomes.size());
			writer.println("Coverable_Peptides\t" + coverable.size());
			writer.println("Coverage_After_Filter\t" + (allMarkerPeptides.isEmpty()
					? 0.0 : (double) coverable.size() / (double) allMarkerPeptides.size()));
			writer.println("Peptides_With_Lost_Marker_Mapping\t" + lostPeptideCount);
			writer.println("FilteredOut_Peptides\t" + filteredOutList.size());
			writer.println("NoUnderscore_Proteins\t" + noUnderscoreProteins.size());
			writer.println("NoUnderscore_MarkerProteins\t" + noUnderscoreMarkerProteins.size());
			writer.println("NoUnderscore_Marker_Examples\t" + String.join(",", noUnderscoreMarkerProteins));
			writer.println("Lost_Marker_Peptide_Examples\t" + String.join(",", lostExamples));
			writer.println();
			writer.println("FilteredOut_Peptides_List\tGenomeCategories");
			for (String peptide : filteredOutList) {
				String union = formatGenomeCategories(peptideUnionGenomeCategories.get(peptide));
				writer.println(peptide + "\t" + union);
			}
			writer.println();
			writer.println("HardFiltered_Genomes_List");
			ArrayList<String> hardFilteredList = new ArrayList<>(hardFilteredGenomes);
			Collections.sort(hardFilteredList);
			for (String genome : hardFilteredList) {
				writer.println(genome);
			}
		}

		System.out.println("MarkerGenomeSelector edge-case test summary:");
		System.out.println("  Total marker peptides: " + allMarkerPeptides.size());
		System.out.println("  Coverable after hard filter: " + coverable.size());
		System.out.println("  Lost marker mapping peptides: " + lostPeptideCount);
		System.out.println("  Filtered-out peptides (non-coverable): " + filteredOutList.size());
		System.out.println("  No-underscore marker proteins: " + noUnderscoreMarkerProteins.size());
		System.out.println("Test outputs:");
		System.out.println("  Diagnostics TSV: " + diagnosticsTsv.getAbsolutePath());
	}

	private Map<String, Set<String>> copyGenomeCategories(Map<String, Set<String>> input) {
		Map<String, Set<String>> copy = new HashMap<>();
		for (Map.Entry<String, Set<String>> entry : input.entrySet()) {
			copy.put(entry.getKey(), new HashSet<>(entry.getValue()));
		}
		return copy;
	}

	private String formatGenomeCategories(Map<String, Set<String>> genomeCategories) {
		if (genomeCategories == null || genomeCategories.isEmpty()) {
			return "";
		}
		ArrayList<String> genomes = new ArrayList<>(genomeCategories.keySet());
		Collections.sort(genomes);
		ArrayList<String> parts = new ArrayList<>();
		for (String genome : genomes) {
			ArrayList<String> cats = new ArrayList<>(genomeCategories.get(genome));
			Collections.sort(cats);
			parts.add(genome + ":" + String.join(",", cats));
		}
		return String.join("|", parts);
	}

	private void runMarkerGenomeSelector(File markerFasta, File peptideTsv, File taxaTsv,
			File evidenceTsv, File rankedTsv, File selectedTsv) {
		try {
			Class<?> selectorClass = Class.forName("uk.ac.quadram.metapilot.task.mag.MarkerGenomeSelector");
			Object selector = selectorClass.getDeclaredConstructor().newInstance();
			java.lang.reflect.Method method = selectorClass.getMethod(
					"peptideTestWithStrategy",
					String.class,
					File.class,
					File.class,
					String.class,
					String.class,
					String.class,
					double.class,
					double.class,
					int.class);
			method.invoke(selector,
					peptideTsv.getAbsolutePath(),
					markerFasta,
					taxaTsv,
					evidenceTsv.getAbsolutePath(),
					rankedTsv.getAbsolutePath(),
					selectedTsv.getAbsolutePath(),
					1.0,
					0.00001,
					200);
		} catch (Exception e) {
			System.out.println("Failed to run MarkerGenomeSelector: " + e.getMessage());
		}
	}
	*/
	public void getGenomeCountFromFasta(String in) {
		HashSet<String> genomeSet = new HashSet<String>();
		try (BufferedReader reader = new BufferedReader(new FileReader(in))) {
			String line = null;
			while ((line = reader.readLine()) != null) {
				if (line.startsWith(">")) {
					String genomeName = line.substring(1, line.indexOf("_"));
					genomeSet.add(genomeName);
				}
			}
			reader.close();
		} catch (IOException e) {
		}
		System.out.println(genomeSet.size());
	}
	
	public static void main(String[] args) throws IOException {
		DatasetMarkerAnalysis analyzer = new DatasetMarkerAnalysis();
		
		// Example usage
//		File peptideFile = new File(
//				"C:\\Users\\daz25hak\\Documents\\papers\\MetaPilot_2025\\reference\\metapep\\15projects_final_peptides"
//						+ "\\PXD001573_final_peptides.tsv");
//				File peptideFile = new File("Z:\\PRIDE_dataset\\PXD056930_300drug\\search"
		
		File peptideFile = new File(
				"C:\\Users\\daz25hak\\Documents\\papers\\MetaPilot_2025\\reference\\metapep\\15projects_final_peptides");
		File faaFolder = new File("C:\\Users\\daz25hak\\MetaProAI\\database\\MAG\\human-gut\\v2.0.2\\original_db");
		File outputFile = new File("C:\\Users\\daz25hak\\MetaPilot\\MetaPilot\\test\\analysis");
		
		// Optional: If protein names don't contain genome IDs, provide predicted folder
		// File predictedFolder = new File("path/to/predicted/folder");
		File predictedFolder = new File("C:\\Users\\daz25hak\\MetaProAI\\database\\MAG\\human-gut\\v2.0.2\\predicted");  // Set to null to use protein ID extraction method

		// Razor parameters
		double razorThreshold = 0.99; // 99% peptide coverage
		double stepThreshold = 0.00001;   // 0.001% increment threshold
		int maxGenomes = 500;           // Maximum genomes to select
		/*
		File[] files = peptideFile.listFiles();
		for(int i=0;i<files.length;i++) {
			String name = files[i].getName();
			String project = name;
			if(name.endsWith("tsv")) {
				project = name.split("_")[0];
			} else if(name.endsWith("txt")) {
				project = name.split("_")[1];
			}
			
			analyzer.analyzeDataset(files[i], faaFolder, predictedFolder, outputFile, project,
					razorThreshold, stepThreshold, maxGenomes);
		}
		*/
		
//		analyzer.analyzeDataset(new File("Z:\\PRIDE_dataset\\PXD004039\\MetaPilot\\sage_MAG_search"
//				+ "\\final_peptides.tsv"), faaFolder, predictedFolder, outputFile, "PXD004039_new",
//				razorThreshold, stepThreshold, maxGenomes);
		
//		analyzer.analyzeDatasetFromProteinInfo(
//				new File("Z:\\PRIDE_dataset\\PXD006118_3mix_type\\Run1\\refdb" + "\\lfq.tsv"),
//				new File("Z:\\PRIDE_dataset\\PXD006118_3mix_type\\Mock_Comm_RefDB_V3_info.txt"), outputFile, "3mix_run1",
//				1.00, stepThreshold, maxGenomes);
		
//		analyzer.getGenomeCount("C:\\Users\\daz25hak\\MetaPilot\\MetaPilot\\test\\analysis\\Figure1D");
		
//		analyzer.extractProteinInfoFromFasta(new File("Z:\\PRIDE_dataset\\PXD006118_3mix_type\\Mock_Comm_RefDB_V3.fasta"), 
//				new File("Z:\\PRIDE_dataset\\PXD006118_3mix_type\\Mock_Comm_RefDB_V3_info.txt"));
		
//		analyzer.calculateCorrelationFromPaperData(new File("C:\\Users\\daz25hak\\MetaPilot\\MetaPilot\\test\\analysis"
//				+ "\\paper_data.txt"), new File("C:\\Users\\daz25hak\\MetaPilot\\MetaPilot\\test\\analysis\\Figure2"
//						+ "\\paper_data_correlation.txt"));
		
//		analyzer.calculateTaxaAbundanceFromPeptides(
//				new File("C:\\Users\\daz25hak\\MetaPilot\\MetaPilot\\test\\analysis\\Run1_lfq.tsv"),
//				new File("C:\\Users\\daz25hak\\MetaPilot\\MetaPilot\\test\\analysis\\paper_data.txt"),
//				new File("C:\\Users\\daz25hak\\MetaPilot\\MetaPilot\\test\\analysis\\Figure2"
//						+ "\\taxa_abundance_correlation.txt"));
		
//		analyzer.calculateTaxaAbundanceFromPeptides(
//				new File("Z:\\PRIDE_dataset\\PXD006118_3mix_type\\Run1\\refdb" + "\\run1_lfq.tsv"),
//				new File("C:\\Users\\daz25hak\\Documents\\papers\\MetaPilot_2025\\reference\\3mix_type"
//						+ "\\paper_data.txt"),
//				new File("C:\\Users\\daz25hak\\MetaPilot\\MetaPilot\\test\\analysis\\Figure2\\allpep_correlation.tsv"));
		
//		analyzer.extractMarkerProteinsFromFasta(new File("Z:\\PRIDE_dataset\\PXD006118_3mix_type"
//				+ "\\Mock_Comm_RefDB_V3.fasta"), new File("Z:\\PRIDE_dataset\\PXD006118_3mix_type\\Mock_Comm_RefDB_V3_marker.fasta"));
		
//		analyzer.exportPeptideFasta("C:\\Users\\daz25hak\\MetaProAI\\database\\MAG\\human-gut\\v2.0.2\\4markers_deepdetect.tsv", 
//				"C:\\Users\\daz25hak\\MetaProAI\\database\\MAG\\human-gut\\v2.0.2\\4markers_highpep.fasta", 0.8);
		
//		analyzer.getSimpleTaxaInfo("C:\\Users\\daz25hak\\MetaProAI\\database\\MAG\\human-gut\\v2.0.2\\genomes-all_metadata.tsv", 
//				"C:\\Users\\daz25hak\\MetaProAI\\database\\MAG\\human-gut\\v2.0.2\\genomes-taxa.tsv");
		
//		analyzer.evaluateGenomeSelection(new File("C:\\Users\\daz25hak\\MetaPilot\\MetaPilot\\test\\analysis\\Figure2"
//				+ "\\genomes.txt"), new File("Z:\\PRIDE_dataset\\PXD006118_3mix_type\\Run1\\MetaPilot\\sage_MAG_search_marker4\\Run1_C1_2000ng\\hap"
//						+ "\\hap_genomes.tsv"));
		
//		analyzer.evaluateGenomeSelection(new File("C:\\Users\\daz25hak\\MetaPilot\\MetaPilot\\test\\analysis\\Figure2"
//				+ "\\genomes.txt"), new File("Z:\\PRIDE_dataset\\PXD006118_3mix_type\\Run1\\MetaPilot\\sage_MAG_search\\Run1_C1_2000ng\\hap"
//						+ "\\hap_genomes.tsv"));
		
//		analyzer.testMarkerGenomeSelectorEdgeCases();

		File resultFolder = new File("Z:\\PRIDE_dataset\\PXD006118_3mix_type\\Run1\\MetaPilot\\sage_MAG_search_marker4\\Run1_C1_2000ng\\hap");
		File markerFasta = new File("C:\\Users\\daz25hak\\MetaProAI\\database\\MAG\\human-gut\\v2.0.2\\4markers.fasta");
 		File peptideTsv = new File(resultFolder, "results.sage.tsv");
 		File taxaTsv = new File("C:\\Users\\daz25hak\\MetaProAI\\database\\MAG\\human-gut\\v2.0.2\\genomes-taxa.tsv");
 		File evidenceTsv = new File(resultFolder, "evidence.tsv");
 		File rankedTsv = new File(resultFolder, "ranked.tsv");
 		File selectedTsv = new File(resultFolder, "selected.tsv");
 		File diagnosticsTsv = new File(resultFolder, "selected_diagnostics.tsv");
		
//		analyzer.peptideTest("Z:\\PRIDE_dataset\\PXD006118_3mix_type\\Run1\\MetaPilot\\sage_MAG_search\\Run1_U1_2000ng\\hap"
//				+ "\\results.sage.tsv", markerFasta);
//		analyzer.peptideTest(
//				"Z:\\PRIDE_dataset\\PXD023217_campi\\MetaPilot\\sage_MAG_search\\F07_Fraction3\\hap"
//						+ "\\results.sage.tsv",
//				markerFasta,
//				"Z:\\PRIDE_dataset\\PXD023217_campi\\MetaPilot\\sage_MAG_search\\F07_Fraction3\\hap"
//						+ "\\genomes.marker4.tsv");
		
// 		analyzer.testMarkerGenomeSelectorEdgeCases(markerFasta, peptideTsv, taxaTsv, evidenceTsv, rankedTsv, selectedTsv, diagnosticsTsv);
 		
 		analyzer.getGenomeCountFromFasta("Z:\\PRIDE_dataset\\PXD006118_3mix_type\\Run1\\MetaPilot\\sage_MAG_search\\combined\\combined.fasta");
 		analyzer.getGenomeCountFromFasta("Z:\\PRIDE_dataset\\PXD006118_3mix_type\\Run1\\MetaPilot\\sage_MAG_search_rib\\combined\\combined.fasta");

	}
}

