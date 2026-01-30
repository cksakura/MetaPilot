package analysis;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;
import java.util.regex.Pattern;

/**
 * Standalone analysis tool for generating peptide distribution statistics.
 * This is separate from the main software for paper publication purposes.
 */
public class PeptideCategoryAnalyzer {
	
	// Regex patterns for the 4-category markers (case-insensitive)
	private static final Pattern ANCHOR_PATTERN = Pattern.compile("(?i)\\bribosomal\\s+protein\\b");
	private static final Pattern HIGH_RES_PATTERN = Pattern.compile("(?i)\\b(gyrase|gyr[ab])\\b|\\brna\\s+polymerase\\b|\\bdna\\s+polymerase\\b|\\brpo[abc]\\b|\\bdna[eq]\\b");
	private static final Pattern ABUNDANCE_PATTERN = Pattern.compile("(?i)\\belongation\\s+factor\\b|\\bef[-\\s]?(tu|g|ts|p|4|1a|2)\\b");
	private static final Pattern STABILITY_PATTERN = Pattern.compile("(?i)(aminoacyl[-\\s]?trna|aa[-\\s]?trna|trna)[-\\s]*[^a-zA-Z0-9]{0,2}\\s*(ligase|synthetase)\\b");
	
	/**
	 * Load protein-to-category mapping from .faa files
	 */
	private Map<String, String> loadMarkerProteinCategories(File folder) throws IOException {
		Map<String, String> map = new HashMap<>();
		File[] files = folder.listFiles();
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
	
	private String categorize(String lineLower) {
		if (ANCHOR_PATTERN.matcher(lineLower).find()) return "Anchor";
		if (HIGH_RES_PATTERN.matcher(lineLower).find()) return "HighRes";
		if (ABUNDANCE_PATTERN.matcher(lineLower).find()) return "Abundance";
		if (STABILITY_PATTERN.matcher(lineLower).find()) return "Stability";
		return null;
	}
	
	/**
	 * List all peptide TSV files in a folder
	 */
	private ArrayList<File> listPeptideFiles(File folder) {
		ArrayList<File> list = new ArrayList<>();
		File[] files = folder.listFiles();
		if (files == null) {
			return list;
		}
		for (File f : files) {
			if (f.isFile() && f.getName().endsWith(".tsv")) {
				list.add(f);
			}
		}
		return list;
	}
	
	/**
	 * Extract genome ID from filename (e.g., "MGYG000000005.predicted.tsv" -> "MGYG000000005")
	 */
	private String genomeIdFromFile(String filename) {
		int dot = filename.indexOf('.');
		if (dot <= 0) {
			return null;
		}
		return filename.substring(0, dot);
	}
	
	/**
	 * Inner class to store category-based statistics for a genome
	 */
	private static class CategoryStats {
		String filename;
		int totalPeptides;
		int anchorCount;
		int highResCount;
		int abundanceCount;
		int stabilityCount;
		int combinedCount;
		
		CategoryStats(String filename, int totalPeptides, int anchor, int highRes, int abundance, int stability) {
			this.filename = filename;
			this.totalPeptides = totalPeptides;
			this.anchorCount = anchor;
			this.highResCount = highRes;
			this.abundanceCount = abundance;
			this.stabilityCount = stability;
			this.combinedCount = anchor + highRes + abundance + stability;
		}
		
		double getAnchorPercentage() { return totalPeptides > 0 ? (anchorCount * 100.0 / totalPeptides) : 0; }
		double getHighResPercentage() { return totalPeptides > 0 ? (highResCount * 100.0 / totalPeptides) : 0; }
		double getAbundancePercentage() { return totalPeptides > 0 ? (abundanceCount * 100.0 / totalPeptides) : 0; }
		double getStabilityPercentage() { return totalPeptides > 0 ? (stabilityCount * 100.0 / totalPeptides) : 0; }
		double getCombinedPercentage() { return totalPeptides > 0 ? (combinedCount * 100.0 / totalPeptides) : 0; }
	}
	
	/**
	 * Get peptide distribution statistics for the 4 marker categories.
	 * Similar to protein statistics but for peptides from predicted TSV files.
	 * @param predictedFolder Folder containing *.predicted.tsv files
	 * @param faaFolder Folder containing .faa files (for protein-to-category mapping)
	 * @param outputFile CSV file to save peptide category statistics
	 * @throws IOException
	 */
	public void getPeptideCategoryStatistics(File predictedFolder, File faaFolder, File outputFile) throws IOException {
		// Load protein-to-category mapping
		Map<String, String> proteinToCategory = loadMarkerProteinCategories(faaFolder);
		
		PrintWriter writer = new PrintWriter(outputFile);
		
		ArrayList<File> peptideFiles = listPeptideFiles(predictedFolder);
		if (peptideFiles.isEmpty()) {
			writer.close();
			throw new IOException("No predicted TSV files found in: " + predictedFolder.getPath());
		}
		
		ArrayList<CategoryStats> statsList = new ArrayList<>();
		
		// Total counters
		int totalFiles = 0;
		int totalPeptides = 0;
		int totalAnchor = 0;
		int totalHighRes = 0;
		int totalAbundance = 0;
		int totalStability = 0;
		
		// Collect statistics for each file
		for (File pepFile : peptideFiles) {
			String genomeId = genomeIdFromFile(pepFile.getName());
			if (genomeId == null) {
				continue;
			}
			
			totalFiles++;
			
			// Use sets to count distinct peptides per category
			Set<String> anchorPeptides = new HashSet<>();
			Set<String> highResPeptides = new HashSet<>();
			Set<String> abundancePeptides = new HashSet<>();
			Set<String> stabilityPeptides = new HashSet<>();
			Set<String> allPeptides = new HashSet<>();
			
			try {
				BufferedReader reader = new BufferedReader(new FileReader(pepFile));
				String line = reader.readLine(); // header
				while ((line = reader.readLine()) != null) {
					String[] parts = line.split("\t");
					if (parts.length < 2) {
						continue;
					}
					String proteinId = parts[0].trim();
					String peptide = parts[1].trim();
					
					if (peptide.isEmpty()) {
						continue;
					}
					
					allPeptides.add(peptide);
					
					// Check category for this protein
					String category = proteinToCategory.get(proteinId);
					if (category != null) {
						switch (category) {
							case "Anchor":
								anchorPeptides.add(peptide);
								break;
							case "HighRes":
								highResPeptides.add(peptide);
								break;
							case "Abundance":
								abundancePeptides.add(peptide);
								break;
							case "Stability":
								stabilityPeptides.add(peptide);
								break;
						}
					}
				}
				reader.close();
				
				int peptideCount = allPeptides.size();
				int anchorCount = anchorPeptides.size();
				int highResCount = highResPeptides.size();
				int abundanceCount = abundancePeptides.size();
				int stabilityCount = stabilityPeptides.size();
				
				totalPeptides += peptideCount;
				totalAnchor += anchorCount;
				totalHighRes += highResCount;
				totalAbundance += abundanceCount;
				totalStability += stabilityCount;
				
				statsList.add(new CategoryStats(pepFile.getName(), peptideCount, 
					anchorCount, highResCount, abundanceCount, stabilityCount));
			} catch (IOException e) {
				System.err.println("Error reading " + pepFile.getName() + ": " + e.getMessage());
			}
		}
		
		// Write CSV header
		writer.println("Filename,Total_Peptides,Anchor,HighResolution,Abundance,Stability,Combined," +
			"Anchor%,HighRes%,Abundance%,Stability%,Combined%");
		
		// Write data
		for (CategoryStats stats : statsList) {
			writer.printf("%s,%d,%d,%d,%d,%d,%d,%.2f%%,%.2f%%,%.2f%%,%.2f%%,%.2f%%\n",
				stats.filename, stats.totalPeptides,
				stats.anchorCount, stats.highResCount, stats.abundanceCount, stats.stabilityCount,
				stats.combinedCount,
				stats.getAnchorPercentage(), stats.getHighResPercentage(),
				stats.getAbundancePercentage(), stats.getStabilityPercentage(),
				stats.getCombinedPercentage());
		}
		
		// Calculate overall statistics
		int totalCombined = totalAnchor + totalHighRes + totalAbundance + totalStability;
		double overallAnchorPct = totalPeptides > 0 ? (totalAnchor * 100.0 / totalPeptides) : 0;
		double overallHighResPct = totalPeptides > 0 ? (totalHighRes * 100.0 / totalPeptides) : 0;
		double overallAbundancePct = totalPeptides > 0 ? (totalAbundance * 100.0 / totalPeptides) : 0;
		double overallStabilityPct = totalPeptides > 0 ? (totalStability * 100.0 / totalPeptides) : 0;
		double overallCombinedPct = totalPeptides > 0 ? (totalCombined * 100.0 / totalPeptides) : 0;
		
		writer.println();
		writer.println("=== OVERALL SUMMARY ===");
		writer.printf("Total Files,%d\n", totalFiles);
		writer.printf("Total Peptides,%d\n", totalPeptides);
		writer.println();
		writer.printf("Anchor (Ribosomal),%d,%.2f%%\n", totalAnchor, overallAnchorPct);
		writer.printf("High Resolution (Poly/Gyr),%d,%.2f%%\n", totalHighRes, overallHighResPct);
		writer.printf("Abundance (Elongation),%d,%.2f%%\n", totalAbundance, overallAbundancePct);
		writer.printf("Stability (tRNA),%d,%.2f%%\n", totalStability, overallStabilityPct);
		writer.printf("Combined (All 4),%d,%.2f%%\n", totalCombined, overallCombinedPct);
		
		// Analyze each category separately
		writer.println();
		writer.println("=== CATEGORY ANALYSIS ===");
		
		String[] categoryNames = {"Anchor", "High Resolution", "Abundance", "Stability", "Combined"};
		for (int cat = 0; cat < 5; cat++) {
			final int category = cat;
			ArrayList<CategoryStats> sortedList = new ArrayList<>(statsList);
			
			// Sort by the specific category
			Collections.sort(sortedList, new Comparator<CategoryStats>() {
				@Override
				public int compare(CategoryStats a, CategoryStats b) {
					int countA = getCategoryCount(a, category);
					int countB = getCategoryCount(b, category);
					return Integer.compare(countB, countA); // descending
				}
			});
			
			// Calculate statistics for this category
			double sumCount = 0;
			double sumPct = 0;
			int minCount = Integer.MAX_VALUE;
			int maxCount = 0;
			
			for (CategoryStats stats : sortedList) {
				int count = getCategoryCount(stats, category);
				double pct = getCategoryPercentage(stats, category);
				sumCount += count;
				sumPct += pct;
				if (count < minCount) minCount = count;
				if (count > maxCount) maxCount = count;
			}
			
			double meanCount = statsList.size() > 0 ? sumCount / statsList.size() : 0;
			double meanPct = statsList.size() > 0 ? sumPct / statsList.size() : 0;
			
			writer.println();
			writer.printf("--- %s ---\n", categoryNames[category]);
			writer.printf("Mean Count,%.2f\n", meanCount);
			writer.printf("Mean Percentage,%.2f%%\n", meanPct);
			writer.printf("Range,%d to %d\n", minCount, maxCount);
			writer.println("Top 5 Genomes:");
			for (int i = 0; i < Math.min(5, sortedList.size()); i++) {
				CategoryStats stats = sortedList.get(i);
				int count = getCategoryCount(stats, category);
				double pct = getCategoryPercentage(stats, category);
				writer.printf("  %d. %s: %d peptides (%.2f%%)\n", i+1, stats.filename, count, pct);
			}
		}
		
		writer.close();
		
		// Print summary to console
		System.out.println("Peptide category statistics saved to: " + outputFile.getPath());
		System.out.println("Total files: " + totalFiles);
		System.out.println("Anchor (Ribosomal): " + totalAnchor + " (" + String.format("%.2f%%", overallAnchorPct) + ")");
		System.out.println("High Resolution: " + totalHighRes + " (" + String.format("%.2f%%", overallHighResPct) + ")");
		System.out.println("Abundance (Elongation): " + totalAbundance + " (" + String.format("%.2f%%", overallAbundancePct) + ")");
		System.out.println("Stability (tRNA): " + totalStability + " (" + String.format("%.2f%%", overallStabilityPct) + ")");
		System.out.println("Combined (All 4): " + totalCombined + " (" + String.format("%.2f%%", overallCombinedPct) + ")");
	}
	
	/**
	 * Helper method to get count for a specific category
	 */
	private int getCategoryCount(CategoryStats stats, int category) {
		switch (category) {
			case 0: return stats.anchorCount;
			case 1: return stats.highResCount;
			case 2: return stats.abundanceCount;
			case 3: return stats.stabilityCount;
			case 4: return stats.combinedCount;
			default: return 0;
		}
	}
	
	/**
	 * Helper method to get percentage for a specific category
	 */
	private double getCategoryPercentage(CategoryStats stats, int category) {
		switch (category) {
			case 0: return stats.getAnchorPercentage();
			case 1: return stats.getHighResPercentage();
			case 2: return stats.getAbundancePercentage();
			case 3: return stats.getStabilityPercentage();
			case 4: return stats.getCombinedPercentage();
			default: return 0;
		}
	}
	
	/**
	 * Calculate distinct peptide counts (peptides found in only ONE genome) per category.
	 * Uses batch processing if memory is limited.
	 * @param predictedFolder Folder containing *.predicted.tsv files
	 * @param faaFolder Folder containing .faa files (for protein-to-category mapping)
	 * @param outputFile CSV file to save distinct peptide statistics
	 * @throws IOException
	 */
	public void getDistinctPeptideStatistics(File predictedFolder, File faaFolder, File outputFile) throws IOException {
		// Load protein-to-category mapping
		Map<String, String> proteinToCategory = loadMarkerProteinCategories(faaFolder);
		
		ArrayList<File> peptideFiles = listPeptideFiles(predictedFolder);
		if (peptideFiles.isEmpty()) {
			throw new IOException("No predicted TSV files found in: " + predictedFolder.getPath());
		}
		
		System.out.println("Found " + peptideFiles.size() + " peptide files.");
		System.out.println("Attempting to calculate distinct peptides (peptides in only one genome)...");
		
		calculateDistinctPeptidesBatch(peptideFiles, proteinToCategory, outputFile);
	}
	
	/**
	 * Single-pass calculation (loads all data into memory)
	 */
	private void calculateDistinctPeptidesSinglePass(ArrayList<File> peptideFiles, 
			Map<String, String> proteinToCategory, File outputFile) throws IOException {
		
		// Map: peptide -> set of genomes that contain it
		Map<String, Set<String>> peptideToGenomes = new HashMap<>();
		// Map: peptide -> category (for marker peptides)
		Map<String, String> peptideToCategory = new HashMap<>();
		
		System.out.println("First pass: collecting all peptides...");
		int processed = 0;
		for (File pepFile : peptideFiles) {
			String genomeId = genomeIdFromFile(pepFile.getName());
			if (genomeId == null) {
				continue;
			}
			
			try {
				BufferedReader reader = new BufferedReader(new FileReader(pepFile));
				String line = reader.readLine(); // header
				while ((line = reader.readLine()) != null) {
					String[] parts = line.split("\t");
					if (parts.length < 2) {
						continue;
					}
					String proteinId = parts[0].trim();
					String peptide = parts[1].trim();
					
					if (peptide.isEmpty()) {
						continue;
					}
					
					// Track which genomes have this peptide
					peptideToGenomes.computeIfAbsent(peptide, k -> new HashSet<>()).add(genomeId);
					
					// Track category for marker peptides
					String category = proteinToCategory.get(proteinId);
					if (category != null && !peptideToCategory.containsKey(peptide)) {
						peptideToCategory.put(peptide, category);
					}
				}
				reader.close();
				
				processed++;
				if (processed % 100 == 0) {
					System.out.println("Processed " + processed + " / " + peptideFiles.size() + " files...");
				}
			} catch (IOException e) {
				System.err.println("Error reading " + pepFile.getName() + ": " + e.getMessage());
			}
		}
		
		System.out.println("Second pass: counting distinct peptides per genome...");
		// Now count distinct peptides per genome per category
		Map<String, Map<String, Integer>> genomeCategoryCounts = new HashMap<>();
		Map<String, Integer> genomeTotalCounts = new HashMap<>(); // Total distinct peptides (all proteins)
		
		for (Map.Entry<String, Set<String>> entry : peptideToGenomes.entrySet()) {
			String peptide = entry.getKey();
			Set<String> genomes = entry.getValue();
			
			// Only count if peptide appears in exactly ONE genome
			if (genomes.size() == 1) {
				String genomeId = genomes.iterator().next();
				
				// Count ALL distinct peptides (from all proteins)
				genomeTotalCounts.merge(genomeId, 1, Integer::sum);
				
				// Count marker-4 distinct peptides
				String category = peptideToCategory.get(peptide);
				if (category != null) {
					genomeCategoryCounts.computeIfAbsent(genomeId, k -> new HashMap<>())
						.merge(category, 1, Integer::sum);
				}
			}
		}
		
		// Write results
		writeDistinctPeptideResults(genomeCategoryCounts, genomeTotalCounts, outputFile);
	}
	
	/**
	 * Batch processing by peptide prefix (A-Z) - writes temp files to save memory
	 */
	private void calculateDistinctPeptidesBatch(ArrayList<File> peptideFiles, 
			Map<String, String> proteinToCategory, File outputFile) throws IOException {
		
		// Process in batches by first character of peptide (amino acids are A-Z)
		String[] prefixes = {"A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", 
			"N", "O", "P", "Q", "R", "S", "T", "U", "V", "W", "X", "Y", "Z"};
		
		// Create temp directory
		File tempDir = new File(outputFile.getParent(), "temp_distinct_" + System.currentTimeMillis());
		tempDir.mkdirs();
		
		ArrayList<File> tempFiles = new ArrayList<>();
		
		try {
			// Process each prefix batch and write to temp file
			for (String prefix : prefixes) {
				System.out.println("Processing peptides starting with '" + prefix + "'...");
				System.gc(); // Suggest garbage collection
				
				// Map: peptide -> set of genomes (only for peptides starting with prefix)
				Map<String, Set<String>> peptideToGenomes = new HashMap<>();
				Map<String, String> peptideToCategory = new HashMap<>();
				
				// First pass: collect peptides starting with this prefix
				for (File pepFile : peptideFiles) {
					String genomeId = genomeIdFromFile(pepFile.getName());
					if (genomeId == null) {
						continue;
					}
					
					try {
						BufferedReader reader = new BufferedReader(new FileReader(pepFile));
						String line = reader.readLine(); // header
						while ((line = reader.readLine()) != null) {
							String[] parts = line.split("\t");
							if (parts.length < 2) {
								continue;
							}
							String proteinId = parts[0].trim();
							String peptide = parts[1].trim();
							
							if (peptide.isEmpty()) {
								continue;
							}
							
							// Case-insensitive prefix matching
							String peptideUpper = peptide.toUpperCase();
							String prefixUpper = prefix.toUpperCase();
							if (!peptideUpper.startsWith(prefixUpper)) {
								continue;
							}
							
							peptideToGenomes.computeIfAbsent(peptide, k -> new HashSet<>()).add(genomeId);
							
							String category = proteinToCategory.get(proteinId);
							if (category != null && !peptideToCategory.containsKey(peptide)) {
								peptideToCategory.put(peptide, category);
							}
						}
						reader.close();
					} catch (IOException e) {
						System.err.println("Error reading " + pepFile.getName() + ": " + e.getMessage());
					}
				}
				
				// Count distinct peptides for this prefix
				// Track both total distinct peptides (all proteins) and marker-4 distinct peptides
				Map<String, Map<String, Integer>> batchCounts = new HashMap<>();
				Map<String, Integer> batchTotalCounts = new HashMap<>(); // genome -> total distinct peptides
				
				for (Map.Entry<String, Set<String>> entry : peptideToGenomes.entrySet()) {
					String peptide = entry.getKey();
					Set<String> genomes = entry.getValue();
					
					if (genomes.size() == 1) {
						String genomeId = genomes.iterator().next();
						
						// Count ALL distinct peptides (from all proteins)
						batchTotalCounts.merge(genomeId, 1, Integer::sum);
						
						// Count marker-4 distinct peptides
						String category = peptideToCategory.get(peptide);
						if (category != null) {
							batchCounts.computeIfAbsent(genomeId, k -> new HashMap<>())
								.merge(category, 1, Integer::sum);
						}
					}
				}
				
				// Write this batch to temp file (include total counts)
				File tempFile = new File(tempDir, "batch_" + prefix + ".tsv");
				File tempTotalFile = new File(tempDir, "batch_" + prefix + "_total.tsv");
				writeBatchResults(batchCounts, tempFile);
				writeBatchTotalResults(batchTotalCounts, tempTotalFile);
				tempFiles.add(tempFile);
				tempFiles.add(tempTotalFile);
				
				// Clear memory for this batch
				peptideToGenomes.clear();
				peptideToCategory.clear();
				batchCounts.clear();
				System.gc();
			}
			
			// Combine all temp files into final result
			System.out.println("Combining results from " + tempFiles.size() + " batches...");
			combineBatchResults(tempFiles, outputFile);
			
		} finally {
			// Clean up temp files
			System.out.println("Cleaning up temporary files...");
			for (File tempFile : tempFiles) {
				if (tempFile.exists()) {
					tempFile.delete();
				}
			}
			if (tempDir.exists()) {
				tempDir.delete();
			}
		}
	}
	
	/**
	 * Write batch results to temp file (TSV format: GenomeId\tCategory\tCount)
	 */
	private void writeBatchResults(Map<String, Map<String, Integer>> batchCounts, File tempFile) throws IOException {
		PrintWriter writer = new PrintWriter(tempFile);
		for (Map.Entry<String, Map<String, Integer>> entry : batchCounts.entrySet()) {
			String genomeId = entry.getKey();
			Map<String, Integer> counts = entry.getValue();
			for (Map.Entry<String, Integer> catEntry : counts.entrySet()) {
				writer.println(genomeId + "\t" + catEntry.getKey() + "\t" + catEntry.getValue());
			}
		}
		writer.close();
	}
	
	/**
	 * Write batch total counts to temp file (TSV format: GenomeId\tTotal)
	 */
	private void writeBatchTotalResults(Map<String, Integer> batchTotalCounts, File tempFile) throws IOException {
		PrintWriter writer = new PrintWriter(tempFile);
		for (Map.Entry<String, Integer> entry : batchTotalCounts.entrySet()) {
			writer.println(entry.getKey() + "\t" + entry.getValue());
		}
		writer.close();
	}
	
	/**
	 * Combine all batch temp files into final result
	 */
	private void combineBatchResults(ArrayList<File> tempFiles, File outputFile) throws IOException {
		// Aggregate counts: genome -> category -> total count
		Map<String, Map<String, Integer>> allCounts = new HashMap<>();
		// Aggregate total distinct peptides: genome -> total count
		Map<String, Integer> allTotalCounts = new HashMap<>();
		
		// Read each temp file and aggregate
		for (File tempFile : tempFiles) {
			BufferedReader reader = new BufferedReader(new FileReader(tempFile));
			String line;
			if (tempFile.getName().contains("_total.tsv")) {
				// This is a total count file
				while ((line = reader.readLine()) != null) {
					String[] parts = line.split("\t");
					if (parts.length != 2) {
						continue;
					}
					String genomeId = parts[0];
					int count = Integer.parseInt(parts[1]);
					allTotalCounts.merge(genomeId, count, Integer::sum);
				}
			} else {
				// This is a category count file
				while ((line = reader.readLine()) != null) {
					String[] parts = line.split("\t");
					if (parts.length != 3) {
						continue;
					}
					String genomeId = parts[0];
					String category = parts[1];
					int count = Integer.parseInt(parts[2]);
					
					allCounts.computeIfAbsent(genomeId, k -> new HashMap<>())
						.merge(category, count, Integer::sum);
				}
			}
			reader.close();
		}
		
		// Write final results
		writeDistinctPeptideResults(allCounts, allTotalCounts, outputFile);
	}
	
	/**
	 * Write distinct peptide statistics to CSV file
	 */
	private void writeDistinctPeptideResults(Map<String, Map<String, Integer>> genomeCategoryCounts, 
			Map<String, Integer> genomeTotalCounts, File outputFile) throws IOException {
		
		PrintWriter writer = new PrintWriter(outputFile);
		
		// Write header
		writer.println("Filename,Total_Distinct_Peptides,Anchor,HighResolution,Abundance,Stability,Combined," +
			"Anchor%,HighRes%,Abundance%,Stability%,Combined%");
		
		// Get all unique genomes (from both maps)
		Set<String> allGenomes = new HashSet<>();
		allGenomes.addAll(genomeCategoryCounts.keySet());
		allGenomes.addAll(genomeTotalCounts.keySet());
		
		// Write data for each genome
		for (String genomeId : allGenomes) {
			Map<String, Integer> counts = genomeCategoryCounts.getOrDefault(genomeId, new HashMap<>());
			int total = genomeTotalCounts.getOrDefault(genomeId, 0); // Real total distinct peptides (all proteins)
			
			int anchor = counts.getOrDefault("Anchor", 0);
			int highRes = counts.getOrDefault("HighRes", 0);
			int abundance = counts.getOrDefault("Abundance", 0);
			int stability = counts.getOrDefault("Stability", 0);
			int combined = anchor + highRes + abundance + stability;
			
			// Percentages are relative to total distinct peptides (all proteins)
			double anchorPct = total > 0 ? (anchor * 100.0 / total) : 0;
			double highResPct = total > 0 ? (highRes * 100.0 / total) : 0;
			double abundancePct = total > 0 ? (abundance * 100.0 / total) : 0;
			double stabilityPct = total > 0 ? (stability * 100.0 / total) : 0;
			double combinedPct = total > 0 ? (combined * 100.0 / total) : 0;
			
			writer.printf("%s,%d,%d,%d,%d,%d,%d,%.2f%%,%.2f%%,%.2f%%,%.2f%%,%.2f%%\n",
				genomeId, total, anchor, highRes, abundance, stability, combined,
				anchorPct, highResPct, abundancePct, stabilityPct, combinedPct);
		}
		
		writer.close();
		System.out.println("Distinct peptide statistics saved to: " + outputFile.getPath());
	}
	
	public static void main(String[] args) throws IOException {
		PeptideCategoryAnalyzer analyzer = new PeptideCategoryAnalyzer();
		
		String dbFolder = "C:\\Users\\daz25hak\\MetaProAI\\database\\MAG\\human-gut\\v2.0.2\\original_db";
		String predictedFolder = "C:\\Users\\daz25hak\\MetaProAI\\database\\MAG\\human-gut\\v2.0.2\\predicted";
		String outputFolder = "C:\\Users\\daz25hak\\MetaPilot\\MetaPilot\\test\\analysis\\";
		
		// Create output folder if it doesn't exist
		File outputDir = new File(outputFolder);
		if (!outputDir.exists()) {
			outputDir.mkdirs();
		}
		
		// Generate peptide distribution statistics (all peptides)
//		analyzer.getPeptideCategoryStatistics(
//				new File(predictedFolder),
//				new File(dbFolder),
//				new File(outputFolder + "humangut_4category_peptides.csv"));
		
		// Generate distinct peptide statistics (peptides in only one genome)
		System.out.println("\n=== Starting distinct peptide analysis ===");
		analyzer.getDistinctPeptideStatistics(
				new File(predictedFolder),
				new File(dbFolder),
				new File(outputFolder + "humangut_4category_distinct_peptides.csv"));
	}
}

