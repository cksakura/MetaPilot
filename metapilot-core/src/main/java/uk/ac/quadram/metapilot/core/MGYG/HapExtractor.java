package uk.ac.quadram.metapilot.core.MGYG;

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
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.concurrent.atomic.AtomicInteger;

public class HapExtractor {

	// Regex patterns for marker categories (case-insensitive)
	private static final Pattern ANCHOR_PATTERN = Pattern.compile("(?i)\\bribosomal\\s+protein\\b");
	private static final Pattern HIGH_RES_PATTERN = Pattern.compile("(?i)\\b(gyrase|gyr[ab])\\b|\\brna\\s+polymerase\\b|\\bdna\\s+polymerase\\b|\\brpo[abc]\\b|\\bdna[eq]\\b");
	private static final Pattern ABUNDANCE_PATTERN = Pattern.compile("(?i)\\belongation\\s+factor\\b|\\bef[-\\s]?(tu|g|ts|p|4|1a|2)\\b");
	private static final Pattern STABILITY_PATTERN = Pattern.compile("(?i)(aminoacyl[-\\s]?trna|aa[-\\s]?trna|trna)[-\\s]*[^a-zA-Z0-9]{0,2}\\s*(ligase|synthetase)\\b");

	public HapExtractor() {
	}
	
	/**
	 * Extract 2-markers (Anchor = ribosomal protein, Abundance = elongation factor) to a FASTA file.
	 * Uses {@link #ANCHOR_PATTERN} and {@link #ABUNDANCE_PATTERN} only. If a genome has no matching
	 * proteins and fewer than 1000 total proteins, the whole genome is written (fallback).
	 *
	 * @param folder folder containing .faa/.fasta genome files
	 * @param output output FASTA file (e.g. rib_elon.fasta)
	 */
	public void extractTwoMarkers(File folder, File output) throws IOException {
		PrintWriter writer = new PrintWriter(output);
		File[] files = folder.listFiles();
		for (int i = 0; i < files.length; i++) {
			if (files[i].getName().endsWith(".faa") || files[i].getName().endsWith(".fasta")) {
				int[] proCounts = new int[2];
				BufferedReader reader = new BufferedReader(new FileReader(files[i]));
				String lineString = null;
				boolean write = false;
				while ((lineString = reader.readLine()) != null) {
					if (lineString.startsWith(">")) {
						proCounts[0]++;
						String lineLower = lineString.toLowerCase();
						boolean find = isAnchor(lineLower) || isAbundance(lineLower);

						if (find) {
							proCounts[1]++;
							write = true;
							writer.println(lineString);
						} else {
							write = false;
						}

					} else {
						if (write) {
							writer.println(lineString);
						}
					}
				}
				reader.close();

				if (proCounts[1] == 0 && proCounts[0] < 1000) {
					reader = new BufferedReader(new FileReader(files[i]));
					while ((lineString = reader.readLine()) != null) {
						writer.println(lineString);
					}
					reader.close();
				}
			}
		}
		writer.close();
	}
	
	/**
	 * Inner class to store genome statistics
	 */
	private static class GenomeStats {
		String filename;
		int totalProteins;
		int matchingProteins;
		double matchPercentage;
		
		GenomeStats(String filename, int totalProteins, int matchingProteins) {
			this.filename = filename;
			this.totalProteins = totalProteins;
			this.matchingProteins = matchingProteins;
			this.matchPercentage = totalProteins > 0 ? (matchingProteins * 100.0 / totalProteins) : 0.0;
		}
	}
	
	/**
	 * Get statistics on matching proteins in each .faa file without extracting sequences
	 * @param folder Folder containing .faa files
	 * @param outputFile CSV file to save statistics (filename, total_proteins, matching_proteins)
	 * @throws IOException
	 */
	public void getStatistics(File folder, File outputFile) throws IOException {
		PrintWriter writer = new PrintWriter(outputFile);
		
		File[] files = folder.listFiles();
		if (files == null) {
			writer.close();
			throw new IOException("Cannot read folder: " + folder.getPath());
		}
		
		ArrayList<GenomeStats> statsList = new ArrayList<>();
		int totalProteins = 0;
		int totalMatching = 0;
		
		// Collect statistics for each file
		for (int i = 0; i < files.length; i++) {
			if (files[i].isFile() && (files[i].getName().endsWith(".faa") || files[i].getName().endsWith(".fasta"))) {
				int proteinCount = 0;
				int matchingCount = 0;
				
				BufferedReader reader = new BufferedReader(new FileReader(files[i]));
				String lineString = null;
				
				while ((lineString = reader.readLine()) != null) {
					if (lineString.startsWith(">")) {
						proteinCount++;
						String lineLower = lineString.toLowerCase();
						if (isAnchor(lineLower) || isAbundance(lineLower)) {
							matchingCount++;
						}
					}
				}
				reader.close();
				
				totalProteins += proteinCount;
				totalMatching += matchingCount;
				
				statsList.add(new GenomeStats(files[i].getName(), proteinCount, matchingCount));
			}
		}
		
		// Sort by match percentage (descending)
		Collections.sort(statsList, new Comparator<GenomeStats>() {
			@Override
			public int compare(GenomeStats a, GenomeStats b) {
				return Double.compare(b.matchPercentage, a.matchPercentage);
			}
		});
		
		// Calculate statistical measures
		double sumPercentages = 0;
		double minPercentage = Double.MAX_VALUE;
		double maxPercentage = 0;
		
		for (GenomeStats stats : statsList) {
			sumPercentages += stats.matchPercentage;
			if (stats.matchPercentage < minPercentage) minPercentage = stats.matchPercentage;
			if (stats.matchPercentage > maxPercentage) maxPercentage = stats.matchPercentage;
		}
		
		double meanPercentage = statsList.size() > 0 ? sumPercentages / statsList.size() : 0;
		
		// Calculate median
		double medianPercentage = 0;
		if (statsList.size() > 0) {
			int mid = statsList.size() / 2;
			if (statsList.size() % 2 == 0) {
				medianPercentage = (statsList.get(mid - 1).matchPercentage + statsList.get(mid).matchPercentage) / 2.0;
			} else {
				medianPercentage = statsList.get(mid).matchPercentage;
			}
		}
		
		// Calculate standard deviation
		double sumSquaredDiff = 0;
		for (GenomeStats stats : statsList) {
			double diff = stats.matchPercentage - meanPercentage;
			sumSquaredDiff += diff * diff;
		}
		double stdDeviation = statsList.size() > 0 ? Math.sqrt(sumSquaredDiff / statsList.size()) : 0;
		
		// Write CSV header
		writer.println("Filename,Total_Proteins,Matching_Proteins,Match_Percentage");
		
		// Write sorted data
		for (GenomeStats stats : statsList) {
			writer.printf("%s,%d,%d,%.2f%%\n", 
				stats.filename, stats.totalProteins, stats.matchingProteins, stats.matchPercentage);
		}
		
		// Write summary statistics
		double overallPercentage = totalProteins > 0 ? (totalMatching * 100.0 / totalProteins) : 0.0;
		writer.println();
		writer.println("=== SUMMARY STATISTICS ===");
		writer.printf("Total Files,%d\n", statsList.size());
		writer.printf("Total Proteins,%d\n", totalProteins);
		writer.printf("Total Matching,%d\n", totalMatching);
		writer.printf("Overall Match Rate,%.2f%%\n", overallPercentage);
		writer.println();
		writer.println("=== DISTRIBUTION ANALYSIS ===");
		writer.printf("Mean Match Rate,%.2f%%\n", meanPercentage);
		writer.printf("Median Match Rate,%.2f%%\n", medianPercentage);
		writer.printf("Std Deviation,%.2f%%\n", stdDeviation);
		writer.printf("Min Match Rate,%.2f%%\n", minPercentage);
		writer.printf("Max Match Rate,%.2f%%\n", maxPercentage);
		writer.printf("Range,%.2f%%\n", maxPercentage - minPercentage);
		
		// Identify high and low preference genomes by RATE (top/bottom 10)
		writer.println();
		writer.println("=== TOP 10 GENOMES BY MATCH RATE (Highest Percentage) ===");
		int topCount = Math.min(10, statsList.size());
		for (int i = 0; i < topCount; i++) {
			GenomeStats stats = statsList.get(i);
			writer.printf("%d. %s,%.2f%%,%d/%d\n", 
				i + 1, stats.filename, stats.matchPercentage, 
				stats.matchingProteins, stats.totalProteins);
		}
		
		writer.println();
		writer.println("=== BOTTOM 10 GENOMES BY MATCH RATE (Lowest Percentage) ===");
		int bottomCount = Math.min(10, statsList.size());
		for (int i = 0; i < bottomCount; i++) {
			GenomeStats stats = statsList.get(statsList.size() - bottomCount + i);
			writer.printf("%d. %s,%.2f%%,%d/%d\n", 
				i + 1, stats.filename, stats.matchPercentage, 
				stats.matchingProteins, stats.totalProteins);
		}
		
		// Sort by absolute count for count-based analysis
		ArrayList<GenomeStats> statsByCount = new ArrayList<>(statsList);
		Collections.sort(statsByCount, new Comparator<GenomeStats>() {
			@Override
			public int compare(GenomeStats a, GenomeStats b) {
				return Integer.compare(b.matchingProteins, a.matchingProteins);
			}
		});
		
		// Calculate count-based statistics
		double sumCounts = 0;
		int minCount = Integer.MAX_VALUE;
		int maxCount = 0;
		
		for (GenomeStats stats : statsByCount) {
			sumCounts += stats.matchingProteins;
			if (stats.matchingProteins < minCount) minCount = stats.matchingProteins;
			if (stats.matchingProteins > maxCount) maxCount = stats.matchingProteins;
		}
		
		double meanCount = statsByCount.size() > 0 ? sumCounts / statsByCount.size() : 0;
		
		// Calculate median count
		double medianCount = 0;
		if (statsByCount.size() > 0) {
			int mid = statsByCount.size() / 2;
			if (statsByCount.size() % 2 == 0) {
				medianCount = (statsByCount.get(mid - 1).matchingProteins + statsByCount.get(mid).matchingProteins) / 2.0;
			} else {
				medianCount = statsByCount.get(mid).matchingProteins;
			}
		}
		
		// Calculate standard deviation for counts
		double sumSquaredDiffCount = 0;
		for (GenomeStats stats : statsByCount) {
			double diff = stats.matchingProteins - meanCount;
			sumSquaredDiffCount += diff * diff;
		}
		double stdDeviationCount = statsByCount.size() > 0 ? Math.sqrt(sumSquaredDiffCount / statsByCount.size()) : 0;
		
		// Write count-based analysis
		writer.println();
		writer.println("=== ABSOLUTE COUNT ANALYSIS ===");
		writer.printf("Mean Matching Count,%.2f\n", meanCount);
		writer.printf("Median Matching Count,%.2f\n", medianCount);
		writer.printf("Std Deviation,%.2f\n", stdDeviationCount);
		writer.printf("Min Matching Count,%d\n", minCount);
		writer.printf("Max Matching Count,%d\n", maxCount);
		writer.printf("Range,%d\n", maxCount - minCount);
		
		writer.println();
		writer.println("=== TOP 10 GENOMES BY COUNT (Most Matching Proteins) ===");
		int topCountByNum = Math.min(10, statsByCount.size());
		for (int i = 0; i < topCountByNum; i++) {
			GenomeStats stats = statsByCount.get(i);
			writer.printf("%d. %s,%d proteins,%.2f%%,total=%d\n", 
				i + 1, stats.filename, stats.matchingProteins, 
				stats.matchPercentage, stats.totalProteins);
		}
		
		writer.println();
		writer.println("=== BOTTOM 10 GENOMES BY COUNT (Fewest Matching Proteins) ===");
		int bottomCountByNum = Math.min(10, statsByCount.size());
		for (int i = 0; i < bottomCountByNum; i++) {
			GenomeStats stats = statsByCount.get(statsByCount.size() - bottomCountByNum + i);
			writer.printf("%d. %s,%d proteins,%.2f%%,total=%d\n", 
				i + 1, stats.filename, stats.matchingProteins, 
				stats.matchPercentage, stats.totalProteins);
		}
		
		writer.close();
		
		// Print summary to console
		System.out.println("Statistics saved to: " + outputFile.getPath());
		System.out.println("Total files processed: " + statsList.size());
		System.out.println("Total proteins: " + totalProteins);
		System.out.println("Total matching proteins: " + totalMatching);
		System.out.printf("Overall match rate: %.2f%%\n", overallPercentage);
		System.out.printf("Mean match rate: %.2f%% (±%.2f%%)\n", meanPercentage, stdDeviation);
		System.out.printf("Range: %.2f%% - %.2f%%\n", minPercentage, maxPercentage);
	}
	
	/**
	 * Inner class to store category-based statistics for a genome
	 */
	private static class CategoryStats {
		String filename;
		int totalProteins;
		int anchorCount;        // ribosomal protein
		int highResCount;       // polymerase, gyrase
		int abundanceCount;     // elongation factor
		int stabilityCount;     // tRNA ligase/synthetase
		int combinedCount;      // all 4 categories
		
		CategoryStats(String filename, int totalProteins, int anchor, int highRes, int abundance, int stability) {
			this.filename = filename;
			this.totalProteins = totalProteins;
			this.anchorCount = anchor;
			this.highResCount = highRes;
			this.abundanceCount = abundance;
			this.stabilityCount = stability;
			this.combinedCount = anchor + highRes + abundance + stability;
		}
		
		double getAnchorPercentage() { return totalProteins > 0 ? (anchorCount * 100.0 / totalProteins) : 0; }
		double getHighResPercentage() { return totalProteins > 0 ? (highResCount * 100.0 / totalProteins) : 0; }
		double getAbundancePercentage() { return totalProteins > 0 ? (abundanceCount * 100.0 / totalProteins) : 0; }
		double getStabilityPercentage() { return totalProteins > 0 ? (stabilityCount * 100.0 / totalProteins) : 0; }
		double getCombinedPercentage() { return totalProteins > 0 ? (combinedCount * 100.0 / totalProteins) : 0; }
	}
	
	/**
	 * Get statistics for the 4 marker categories (matching Python analysis)
	 * @param folder Folder containing .faa files
	 * @param outputFile CSV file to save category statistics
	 * @throws IOException
	 */
	public void getCategoryStatistics(File folder, File outputFile) throws IOException {
		PrintWriter writer = new PrintWriter(outputFile);
		
		File[] files = folder.listFiles();
		if (files == null) {
			writer.close();
			throw new IOException("Cannot read folder: " + folder.getPath());
		}
		
		ArrayList<CategoryStats> statsList = new ArrayList<>();
		
		// Total counters
		int totalFiles = 0;
		int totalProteins = 0;
		int totalAnchor = 0;
		int totalHighRes = 0;
		int totalAbundance = 0;
		int totalStability = 0;
		
		// Collect statistics for each file
		for (int i = 0; i < files.length; i++) {
			if (files[i].isFile() && (files[i].getName().endsWith(".faa") || files[i].getName().endsWith(".fasta"))) {
				totalFiles++;
				
				int proteinCount = 0;
				int anchorCount = 0;
				int highResCount = 0;
				int abundanceCount = 0;
				int stabilityCount = 0;
				
				BufferedReader reader = new BufferedReader(new FileReader(files[i]));
				String lineString = null;
				
				while ((lineString = reader.readLine()) != null) {
					if (lineString.startsWith(">")) {
						proteinCount++;
						String lineLower = lineString.toLowerCase();
						
						// Check each category (case-insensitive)
						// Anchor: ribosomal protein
						if (isAnchor(lineLower)) {
							anchorCount++;
						}
						
						// High Resolution: polymerase or gyrase
						if (isHighRes(lineLower)) {
							highResCount++;
						}
						
						// Abundance: elongation factor
						if (isAbundance(lineLower)) {
							abundanceCount++;
						}
						
						// Stability: tRNA ligase or synthetase
						if (isStability(lineLower)) {
							stabilityCount++;
						}
					}
				}
				reader.close();
				
				totalProteins += proteinCount;
				totalAnchor += anchorCount;
				totalHighRes += highResCount;
				totalAbundance += abundanceCount;
				totalStability += stabilityCount;
				
				statsList.add(new CategoryStats(files[i].getName(), proteinCount, 
					anchorCount, highResCount, abundanceCount, stabilityCount));
			}
		}
		
		// Write CSV header
		writer.println("Filename,Total_Proteins,Anchor,HighResolution,Abundance,Stability,Combined," +
			"Anchor%,HighRes%,Abundance%,Stability%,Combined%");
		
		// Write data
		for (CategoryStats stats : statsList) {
			writer.printf("%s,%d,%d,%d,%d,%d,%d,%.2f%%,%.2f%%,%.2f%%,%.2f%%,%.2f%%\n",
				stats.filename, stats.totalProteins,
				stats.anchorCount, stats.highResCount, stats.abundanceCount, stats.stabilityCount,
				stats.combinedCount,
				stats.getAnchorPercentage(), stats.getHighResPercentage(),
				stats.getAbundancePercentage(), stats.getStabilityPercentage(),
				stats.getCombinedPercentage());
		}
		
		// Calculate overall statistics
		int totalCombined = totalAnchor + totalHighRes + totalAbundance + totalStability;
		double overallAnchorPct = totalProteins > 0 ? (totalAnchor * 100.0 / totalProteins) : 0;
		double overallHighResPct = totalProteins > 0 ? (totalHighRes * 100.0 / totalProteins) : 0;
		double overallAbundancePct = totalProteins > 0 ? (totalAbundance * 100.0 / totalProteins) : 0;
		double overallStabilityPct = totalProteins > 0 ? (totalStability * 100.0 / totalProteins) : 0;
		double overallCombinedPct = totalProteins > 0 ? (totalCombined * 100.0 / totalProteins) : 0;
		
		writer.println();
		writer.println("=== OVERALL SUMMARY ===");
		writer.printf("Total Files,%d\n", totalFiles);
		writer.printf("Total Proteins,%d\n", totalProteins);
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
				writer.printf("  %d. %s: %d proteins (%.2f%%)\n", i+1, stats.filename, count, pct);
			}
		}
		
		writer.close();
		
		// Print summary to console
		System.out.println("Category statistics saved to: " + outputFile.getPath());
		System.out.println("Total files: " + totalFiles);
		System.out.println("Anchor (Ribosomal): " + totalAnchor + " (" + String.format("%.2f%%", overallAnchorPct) + ")");
		System.out.println("High Resolution: " + totalHighRes + " (" + String.format("%.2f%%", overallHighResPct) + ")");
		System.out.println("Abundance (Elongation): " + totalAbundance + " (" + String.format("%.2f%%", overallAbundancePct) + ")");
		System.out.println("Stability (tRNA): " + totalStability + " (" + String.format("%.2f%%", overallStabilityPct) + ")");
		System.out.println("Combined (All 4): " + totalCombined + " (" + String.format("%.2f%%", overallCombinedPct) + ")");
	}
	
	/**
	 * Extract 4-markers to a combined FASTA file: Anchor (ribosomal protein), High Resolution
	 * (polymerase/gyrase), Abundance (elongation factor), Stability (tRNA ligase/synthetase).
	 * Use this for the 4markers_highpep workflow; use {@link #extractTwoMarkers(File, File)} for
	 * ribosomal + elongation only (rib_elon).
	 *
	 * @param folder folder containing .faa/.fasta genome files
	 * @param output output FASTA file (e.g. 4markers_highpep.fasta)
	 */
	public void extractFourCategories(File folder, File output) throws IOException {
		PrintWriter writer = new PrintWriter(output);
		
		File[] files = folder.listFiles();
		if (files == null) {
			writer.close();
			throw new IOException("Cannot read folder: " + folder.getPath());
		}
		
		int totalFiles = 0;
		int totalProteins = 0;
		int extractedProteins = 0;
		
		for (int i = 0; i < files.length; i++) {
			if (files[i].isFile() && (files[i].getName().endsWith(".faa") || files[i].getName().endsWith(".fasta"))) {
				totalFiles++;
				
				BufferedReader reader = new BufferedReader(new FileReader(files[i]));
				String lineString = null;
				boolean write = false;
				
				while ((lineString = reader.readLine()) != null) {
					if (lineString.startsWith(">")) {
						totalProteins++;
						String lineLower = lineString.toLowerCase();
						
						// Check if matches any of the 4 categories
						boolean match = false;
						
						// Anchor: ribosomal protein
						if (isAnchor(lineLower)) {
							match = true;
						}
						// High Resolution: polymerase or gyrase
						else if (isHighRes(lineLower)) {
							match = true;
						}
						// Abundance: elongation factor
						else if (isAbundance(lineLower)) {
							match = true;
						}
						// Stability: tRNA ligase or synthetase
						else if (isStability(lineLower)) {
							match = true;
						}
						
						if (match) {
							write = true;
							extractedProteins++;
							writer.println(lineString);
						} else {
							write = false;
						}
					} else {
						if (write) {
							writer.println(lineString);
						}
					}
				}
				reader.close();
			}
		}
		
		writer.close();
		
		// Print summary
		System.out.println("4-Category extraction completed!");
		System.out.println("Output file: " + output.getPath());
		System.out.println("Files processed: " + totalFiles);
		System.out.println("Total proteins: " + totalProteins);
		System.out.println("Extracted proteins: " + extractedProteins);
		System.out.printf("Extraction rate: %.2f%%\n", totalProteins > 0 ? (extractedProteins * 100.0 / totalProteins) : 0);
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
	
	private boolean isAnchor(String lineLower) {
		return ANCHOR_PATTERN.matcher(lineLower).find();
	}
	
	private boolean isHighRes(String lineLower) {
		return HIGH_RES_PATTERN.matcher(lineLower).find();
	}
	
	private boolean isAbundance(String lineLower) {
		return ABUNDANCE_PATTERN.matcher(lineLower).find();
	}
	
	private boolean isStability(String lineLower) {
		return STABILITY_PATTERN.matcher(lineLower).find();
	}
	
	/**
	 * Summarize peptide sharing at category level:
	 * For each marker category (Anchor/HighRes/Abundance/Stability), count distinct peptides
	 * observed in 1, 2, 3, ... genomes. Columns are emitted up to the maximum observed share
	 * to avoid long zero-only tails.
	 *
	 * @param faaFolder folder containing .faa/.fasta genomes
	 * @param peptideFolder folder containing peptide TSVs (e.g., *.predicted.tsv)
	 * @param outputFile CSV to write summary
	 */
	public void summarizePeptideSharingByCategory(File faaFolder, File peptideFolder, File outputFile) throws IOException {
		Map<String, String> proteinToCategory = loadMarkerProteinCategories(faaFolder);
		if (proteinToCategory.isEmpty()) {
			throw new IOException("No marker proteins found in: " + faaFolder.getPath());
		}

		ArrayList<File> peptideFiles = listPeptideFiles(peptideFolder);
		if (peptideFiles.isEmpty()) {
			throw new IOException("No peptide TSV files found in: " + peptideFolder.getPath());
		}

		// category -> peptide -> genomes
		Map<String, Map<String, Set<String>>> categoryPeptideGenomes = new HashMap<>();

		for (File pepFile : peptideFiles) {
			String genomeId = genomeIdFromFile(pepFile.getName());
			if (genomeId == null) {
				continue;
			}
			BufferedReader reader = new BufferedReader(new FileReader(pepFile));
			String line = reader.readLine(); // header
			while ((line = reader.readLine()) != null) {
				String[] parts = line.split("\t");
				if (parts.length < 2) {
					continue;
				}
				String proteinId = parts[0];
				String peptide = parts[1];
				String category = proteinToCategory.get(proteinId);
				if (category == null) {
					continue; // not a marker protein
				}
				Map<String, Set<String>> pepMap = categoryPeptideGenomes.computeIfAbsent(category, k -> new HashMap<>());
				pepMap.computeIfAbsent(peptide, k -> new HashSet<>()).add(genomeId);
			}
			reader.close();
		}

		// Compute share counts per category
		Map<String, Map<Integer, Integer>> categoryShareCounts = new HashMap<>();
		Map<String, Integer> categoryTotals = new HashMap<>();
		int globalMaxShare = 0;

		for (Map.Entry<String, Map<String, Set<String>>> catEntry : categoryPeptideGenomes.entrySet()) {
			String category = catEntry.getKey();
			Map<Integer, Integer> counts = categoryShareCounts.computeIfAbsent(category, k -> new HashMap<>());
			int total = 0;
			for (Set<String> genomes : catEntry.getValue().values()) {
				int share = Math.max(1, genomes.size());
				globalMaxShare = Math.max(globalMaxShare, share);
				counts.put(share, counts.getOrDefault(share, 0) + 1);
				total++;
			}
			categoryTotals.put(category, total);
		}

		// Header up to globalMaxShare (avoids excessive zero columns)
		PrintWriter writer = new PrintWriter(outputFile);
		writer.print("Category,TotalDistinctPeptides");
		for (int i = 1; i <= globalMaxShare; i++) {
			writer.print(",Peptides_in_" + i + "_genomes");
		}
		writer.println();

		for (String category : categoryShareCounts.keySet()) {
			Map<Integer, Integer> counts = categoryShareCounts.get(category);
			int total = categoryTotals.getOrDefault(category, 0);
			writer.print(category);
			writer.print("," + total);
			for (int i = 1; i <= globalMaxShare; i++) {
				writer.print("," + counts.getOrDefault(i, 0));
			}
			writer.println();
		}
		writer.close();

		System.out.println("Category-level peptide sharing summary written to: " + outputFile.getPath());
	}

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
				String lower = line.toLowerCase();
				String category = categorize(lower);
				if (category != null) {
					map.put(proteinId, category);
				}
			}
			reader.close();
		}
		return map;
	}

	/**
	 * Summarize, for each 4-category protein, how many of its distinct peptides are observed
	 * across 1, 2, 3, ... genomes based on in-silico digested TSV files.
	 *
	 * @param faaFolder folder containing .faa/.fasta genomes
	 * @param peptideFolder folder containing peptide TSVs (e.g., *.predicted.tsv)
	 * @param outputFile CSV to write summary
	 */
	public void summarizePeptideSharing(File faaFolder, File peptideFolder, File outputFile) throws IOException {
		Map<String, ProteinInfo> proteinMeta = loadMarkerProteins(faaFolder);
		if (proteinMeta.isEmpty()) {
			throw new IOException("No marker proteins found in: " + faaFolder.getPath());
		}

		ArrayList<File> peptideFiles = listPeptideFiles(peptideFolder);
		if (peptideFiles.isEmpty()) {
			throw new IOException("No peptide TSV files found in: " + peptideFolder.getPath());
		}

		// peptide -> genomes containing it (only counting marker proteins)
		Map<String, Set<String>> peptideGenomes = new HashMap<>();
		for (File pepFile : peptideFiles) {
			String genomeId = genomeIdFromFile(pepFile.getName());
			if (genomeId == null) {
				continue;
			}
			BufferedReader reader = new BufferedReader(new FileReader(pepFile));
			String line = reader.readLine(); // header
			while ((line = reader.readLine()) != null) {
				String[] parts = line.split("\t");
				if (parts.length < 2) {
					continue;
				}
				String proteinId = parts[0];
				String peptide = parts[1];
				if (!proteinMeta.containsKey(proteinId)) {
					continue;
				}
				peptideGenomes.computeIfAbsent(peptide, k -> new HashSet<>()).add(genomeId);
			}
			reader.close();
		}

		// protein -> (shareCount -> peptide count)
		Map<String, Map<Integer, Integer>> proteinShareCounts = new HashMap<>();
		Map<String, Set<String>> proteinSeenPeptides = new HashMap<>();
		int maxShare = 0;

		for (File pepFile : peptideFiles) {
			String genomeId = genomeIdFromFile(pepFile.getName());
			if (genomeId == null) {
				continue;
			}
			BufferedReader reader = new BufferedReader(new FileReader(pepFile));
			String line = reader.readLine(); // header
			while ((line = reader.readLine()) != null) {
				String[] parts = line.split("\t");
				if (parts.length < 2) {
					continue;
				}
				String proteinId = parts[0];
				String peptide = parts[1];
				if (!proteinMeta.containsKey(proteinId)) {
					continue;
				}

				// Deduplicate per protein
				Set<String> seen = proteinSeenPeptides.computeIfAbsent(proteinId, k -> new HashSet<>());
				if (!seen.add(peptide)) {
					continue;
				}

				Set<String> genomes = peptideGenomes.get(peptide);
				int share = genomes == null ? 1 : genomes.size();
				maxShare = Math.max(maxShare, share);

				Map<Integer, Integer> shareCounts = proteinShareCounts.computeIfAbsent(proteinId, k -> new HashMap<>());
				shareCounts.put(share, shareCounts.getOrDefault(share, 0) + 1);
			}
			reader.close();
		}

		// Write summary
		PrintWriter writer = new PrintWriter(outputFile);
		writer.print("ProteinId,GenomeId,Category,ProteinName,TotalDistinctPeptides");
		for (int i = 1; i <= maxShare; i++) {
			writer.print(",Peptides_in_" + i + "_genomes");
		}
		writer.println();

		for (Map.Entry<String, ProteinInfo> entry : proteinMeta.entrySet()) {
			String proteinId = entry.getKey();
			ProteinInfo info = entry.getValue();
			Map<Integer, Integer> counts = proteinShareCounts.getOrDefault(proteinId, new HashMap<>());
			int total = 0;
			for (Integer c : counts.values()) {
				total += c;
			}

			writer.print(proteinId);
			writer.print("," + info.genomeId);
			writer.print("," + info.category);
			writer.print(",\"" + info.name.replace("\"", "") + "\"");
			writer.print("," + total);
			for (int i = 1; i <= maxShare; i++) {
				writer.print("," + counts.getOrDefault(i, 0));
			}
			writer.println();
		}
		writer.close();

		System.out.println("Peptide sharing summary written to: " + outputFile.getPath());
	}

	private Map<String, ProteinInfo> loadMarkerProteins(File folder) throws IOException {
		Map<String, ProteinInfo> map = new HashMap<>();
		File[] files = folder.listFiles();
		if (files == null) {
			return map;
		}
		for (File f : files) {
			if (!(f.isFile() && (f.getName().endsWith(".faa") || f.getName().endsWith(".fasta")))) {
				continue;
			}
			String genomeId = genomeIdFromFile(f.getName());
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
				String name = parts.length > 1 ? parts[1] : "";
				String lower = line.toLowerCase();
				String category = categorize(lower);
				if (category != null) {
					map.put(proteinId, new ProteinInfo(proteinId, genomeId, name, category));
				}
			}
			reader.close();
		}
		return map;
	}

	private String categorize(String lineLower) {
		if (isAnchor(lineLower)) return "Anchor";
		if (isHighRes(lineLower)) return "HighRes";
		if (isAbundance(lineLower)) return "Abundance";
		if (isStability(lineLower)) return "Stability";
		return null;
	}

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

	private String genomeIdFromFile(String filename) {
		int dot = filename.indexOf('.');
		if (dot <= 0) {
			return null;
		}
		return filename.substring(0, dot);
	}

	private static class ProteinInfo {
		String id;
		String genomeId;
		String name;
		String category;

		ProteinInfo(String id, String genomeId, String name, String category) {
			this.id = id;
			this.genomeId = genomeId;
			this.name = name;
			this.category = category;
		}
	}
	
	/**
	 * Generate per-genome, per-category protein and peptide counts for the 4 marker categories.
	 * Output: GenomeId, Category, ProteinCount, PeptideCount
	 *
	 * @param faaFolder folder containing .faa/.fasta genomes
	 * @param peptideFolder folder containing peptide TSVs (e.g., *.predicted.tsv)
	 * @param outputFile TSV to write: GenomeId, Category, ProteinCount, PeptideCount
	 */
	public void generateProteinPeptideDistribution(File faaFolder, File peptideFolder, File outputFile) throws IOException {
		// Step 1: Count marker proteins per category per genome from .faa files
		Map<String, Map<String, Integer>> genomeCategoryProteinCount = new HashMap<>();
		
		File[] faaFiles = faaFolder.listFiles();
		if (faaFiles == null) {
			throw new IOException("Cannot read folder: " + faaFolder.getPath());
		}
		
		for (File faaFile : faaFiles) {
			if (!faaFile.isFile() || (!faaFile.getName().endsWith(".faa") && !faaFile.getName().endsWith(".fasta"))) {
				continue;
			}
			
			String genomeId = genomeIdFromFile(faaFile.getName());
			if (genomeId == null) {
				continue;
			}
			
			Map<String, Integer> categoryCounts = new HashMap<>();
			categoryCounts.put("Anchor", 0);
			categoryCounts.put("HighRes", 0);
			categoryCounts.put("Abundance", 0);
			categoryCounts.put("Stability", 0);
			
			BufferedReader reader = new BufferedReader(new FileReader(faaFile));
			String line;
			while ((line = reader.readLine()) != null) {
				if (line.startsWith(">")) {
					String lineLower = line.toLowerCase();
					String category = null;
					
					if (isAnchor(lineLower)) {
						category = "Anchor";
					} else if (isHighRes(lineLower)) {
						category = "HighRes";
					} else if (isAbundance(lineLower)) {
						category = "Abundance";
					} else if (isStability(lineLower)) {
						category = "Stability";
					}
					
					if (category != null) {
						categoryCounts.put(category, categoryCounts.get(category) + 1);
					}
				}
			}
			reader.close();
			
			genomeCategoryProteinCount.put(genomeId, categoryCounts);
		}
		
		// Step 2: Count distinct peptides per category per genome from peptide TSVs
		Map<String, Map<String, Set<String>>> genomeCategoryPeptides = new HashMap<>();
		
		// First, build protein-to-category mapping
		Map<String, String> proteinToCategory = loadMarkerProteinCategories(faaFolder);
		
		ArrayList<File> peptideFiles = listPeptideFiles(peptideFolder);
		for (File pepFile : peptideFiles) {
			String genomeId = genomeIdFromFile(pepFile.getName());
			if (genomeId == null) {
				continue;
			}
			
			Map<String, Set<String>> categoryPeptides = new HashMap<>();
			categoryPeptides.put("Anchor", new HashSet<>());
			categoryPeptides.put("HighRes", new HashSet<>());
			categoryPeptides.put("Abundance", new HashSet<>());
			categoryPeptides.put("Stability", new HashSet<>());
			
			BufferedReader reader = new BufferedReader(new FileReader(pepFile));
			String line = reader.readLine(); // header
			while ((line = reader.readLine()) != null) {
				String[] parts = line.split("\t");
				if (parts.length < 2) {
					continue;
				}
				String proteinId = parts[0].trim();
				String peptide = parts[1].trim();
				
				String category = proteinToCategory.get(proteinId);
				if (category != null && categoryPeptides.containsKey(category)) {
					categoryPeptides.get(category).add(peptide);
				}
			}
			reader.close();
			
			genomeCategoryPeptides.put(genomeId, categoryPeptides);
		}
		
		// Step 3: Write output table
		PrintWriter writer = new PrintWriter(outputFile);
		writer.println("GenomeId\tCategory\tProteinCount\tPeptideCount");
		
		// Get all unique genome IDs
		Set<String> allGenomes = new HashSet<>();
		allGenomes.addAll(genomeCategoryProteinCount.keySet());
		allGenomes.addAll(genomeCategoryPeptides.keySet());
		
		String[] categories = {"Anchor", "HighRes", "Abundance", "Stability"};
		for (String genomeId : allGenomes) {
			for (String category : categories) {
				int proteinCount = genomeCategoryProteinCount.containsKey(genomeId) 
					? genomeCategoryProteinCount.get(genomeId).getOrDefault(category, 0) : 0;
				int peptideCount = genomeCategoryPeptides.containsKey(genomeId)
					? genomeCategoryPeptides.get(genomeId).getOrDefault(category, new HashSet<>()).size() : 0;
				
				writer.printf("%s\t%s\t%d\t%d%n", genomeId, category, proteinCount, peptideCount);
			}
		}
		
		writer.close();
		System.out.println("Protein and peptide distribution written to: " + outputFile.getPath());
	}
	
	/**
	 * Calculate average protein length per category for marker proteins.
	 * Output: Category, AverageLength, MedianLength, MinLength, MaxLength, TotalProteins
	 *
	 * @param faaFolder folder containing .faa/.fasta genomes
	 * @param outputFile TSV to write statistics
	 */
	public void calculateProteinLengthStats(File faaFolder, File outputFile) throws IOException {
		// category -> list of protein lengths
		Map<String, ArrayList<Integer>> categoryLengths = new HashMap<>();
		categoryLengths.put("Anchor", new ArrayList<>());
		categoryLengths.put("HighRes", new ArrayList<>());
		categoryLengths.put("Abundance", new ArrayList<>());
		categoryLengths.put("Stability", new ArrayList<>());
		
		File[] faaFiles = faaFolder.listFiles();
		if (faaFiles == null) {
			throw new IOException("Cannot read folder: " + faaFolder.getPath());
		}
		
		for (File faaFile : faaFiles) {
			if (!faaFile.isFile() || (!faaFile.getName().endsWith(".faa") && !faaFile.getName().endsWith(".fasta"))) {
				continue;
			}
			
			BufferedReader reader = new BufferedReader(new FileReader(faaFile));
			String line;
			String currentCategory = null;
			StringBuilder sequence = new StringBuilder();
			
			while ((line = reader.readLine()) != null) {
				if (line.startsWith(">")) {
					// Process previous protein if any
					if (currentCategory != null && sequence.length() > 0) {
						categoryLengths.get(currentCategory).add(sequence.length());
					}
					
					// Start new protein
					sequence.setLength(0);
					String lineLower = line.toLowerCase();
					
					if (isAnchor(lineLower)) {
						currentCategory = "Anchor";
					} else if (isHighRes(lineLower)) {
						currentCategory = "HighRes";
					} else if (isAbundance(lineLower)) {
						currentCategory = "Abundance";
					} else if (isStability(lineLower)) {
						currentCategory = "Stability";
					} else {
						currentCategory = null;
					}
				} else if (currentCategory != null) {
					// Accumulate sequence
					sequence.append(line.trim());
				}
			}
			
			// Process last protein
			if (currentCategory != null && sequence.length() > 0) {
				categoryLengths.get(currentCategory).add(sequence.length());
			}
			
			reader.close();
		}
		
		// Calculate statistics and write output
		PrintWriter writer = new PrintWriter(outputFile);
		writer.println("Category\tAverageLength\tMedianLength\tMinLength\tMaxLength\tTotalProteins");
		
		String[] categories = {"Anchor", "HighRes", "Abundance", "Stability"};
		for (String category : categories) {
			ArrayList<Integer> lengths = categoryLengths.get(category);
			if (lengths.isEmpty()) {
				writer.printf("%s\t0\t0\t0\t0\t0%n", category);
				continue;
			}
			
			Collections.sort(lengths);
			int total = lengths.size();
			long sum = 0;
			for (int len : lengths) {
				sum += len;
			}
			double avg = (double) sum / total;
			int median = lengths.get(total / 2);
			if (total % 2 == 0 && total > 1) {
				median = (lengths.get(total / 2 - 1) + lengths.get(total / 2)) / 2;
			}
			int min = lengths.get(0);
			int max = lengths.get(total - 1);
			
			writer.printf("%s\t%.2f\t%d\t%d\t%d\t%d%n", category, avg, median, min, max, total);
		}
		
		writer.close();
		System.out.println("Protein length statistics written to: " + outputFile.getPath());
	}
	
	/**
	 * Count distinct peptides per genome per category - peptides that appear ONLY in that genome
	 * (not shared with any other genome). This shows genome-unique peptides.
	 * Output: GenomeId, Category, DistinctPeptideCount
	 *
	 * @param faaFolder folder containing .faa/.fasta genomes
	 * @param peptideFolder folder containing peptide TSVs (e.g., *.predicted.tsv)
	 * @param outputFile TSV to write: GenomeId, Category, DistinctPeptideCount
	 */
	public void generateDistinctPeptideDistribution(File faaFolder, File peptideFolder, File outputFile) throws IOException {
		// Step 1: Build protein-to-category mapping
		Map<String, String> proteinToCategory = loadMarkerProteinCategories(faaFolder);
		
		// Step 2: Collect all peptides per genome per category
		// genome -> category -> set of peptides
		Map<String, Map<String, Set<String>>> genomeCategoryPeptides = new HashMap<>();
		
		// Also track: peptide -> set of genomes (to identify unique peptides)
		Map<String, Set<String>> peptideToGenomes = new HashMap<>();
		
		ArrayList<File> peptideFiles = listPeptideFiles(peptideFolder);
		for (File pepFile : peptideFiles) {
			String genomeId = genomeIdFromFile(pepFile.getName());
			if (genomeId == null) {
				continue;
			}
			
			Map<String, Set<String>> categoryPeptides = new HashMap<>();
			categoryPeptides.put("Anchor", new HashSet<>());
			categoryPeptides.put("HighRes", new HashSet<>());
			categoryPeptides.put("Abundance", new HashSet<>());
			categoryPeptides.put("Stability", new HashSet<>());
			
			BufferedReader reader = new BufferedReader(new FileReader(pepFile));
			String line = reader.readLine(); // header
			while ((line = reader.readLine()) != null) {
				String[] parts = line.split("\t");
				if (parts.length < 2) {
					continue;
				}
				String proteinId = parts[0].trim();
				String peptide = parts[1].trim();
				
				String category = proteinToCategory.get(proteinId);
				if (category != null && categoryPeptides.containsKey(category)) {
					categoryPeptides.get(category).add(peptide);
					
					// Track which genomes this peptide appears in
					peptideToGenomes.putIfAbsent(peptide, new HashSet<>());
					peptideToGenomes.get(peptide).add(genomeId);
				}
			}
			reader.close();
			
			genomeCategoryPeptides.put(genomeId, categoryPeptides);
		}
		
		// Step 3: Count distinct peptides (appear in only one genome) per genome per category
		PrintWriter writer = new PrintWriter(outputFile);
		writer.println("GenomeId\tCategory\tDistinctPeptideCount");
		
		String[] categories = {"Anchor", "HighRes", "Abundance", "Stability"};
		for (String genomeId : genomeCategoryPeptides.keySet()) {
			Map<String, Set<String>> categoryPeptides = genomeCategoryPeptides.get(genomeId);
			
			for (String category : categories) {
				Set<String> peptides = categoryPeptides.getOrDefault(category, new HashSet<>());
				int distinctCount = 0;
				
				// Count peptides that appear ONLY in this genome
				for (String peptide : peptides) {
					Set<String> genomes = peptideToGenomes.get(peptide);
					if (genomes != null && genomes.size() == 1 && genomes.contains(genomeId)) {
						distinctCount++;
					}
				}
				
				writer.printf("%s\t%s\t%d%n", genomeId, category, distinctCount);
			}
		}
		
		writer.close();
		System.out.println("Distinct peptide distribution written to: " + outputFile.getPath());
	}
	
	/**
	 * Calculate protein abundances for 4 markers from real experimental data matrix.
	 * Uses multi-threading to process 4744 genomes efficiently.
	 * 
	 * Two approaches:
	 * 1. All peptides: sum intensities for all matching peptides (even if shared)
	 * 2. Distinct peptides only: only count peptides unique to that genome
	 * 
	 * @param matrixFile TSV with Stripped.Sequence and columns ending with .mzML (intensities)
	 * @param predictedFolder folder containing *.predicted.tsv files (peptide->protein->genome mapping)
	 * @param faaFolder folder containing .faa/.fasta genomes (for protein category mapping)
	 * @param outputFile TSV output: Sample, GenomeId, Category, Abundance_All, Abundance_Distinct, Rank_All, Rank_Distinct
	 * @param numThreads number of threads for parallel processing (default: CPU cores)
	 */
	public void calculateProteinAbundances(File matrixFile, File predictedFolder, File faaFolder, 
			File outputFile, int numThreads) throws IOException {
		if (numThreads <= 0) {
			numThreads = Runtime.getRuntime().availableProcessors();
		}
		
		System.out.println("Step 1: Loading protein-to-category mapping...");
		Map<String, String> proteinToCategory = loadMarkerProteinCategories(faaFolder);
		
		System.out.println("Step 2: Reading experimental data matrix (peptides with intensities only)...");
		// First, read matrix file and store only peptides that have intensities
		BufferedReader reader = new BufferedReader(new FileReader(matrixFile), 8192 * 4);
		String headerLine = reader.readLine();
		if (headerLine == null) {
			reader.close();
			throw new IOException("Empty matrix file");
		}
		
		String[] headers = headerLine.split("\t");
		int seqColIdx = -1;
		ArrayList<Integer> intensityColIndices = new ArrayList<>();
		
		for (int i = 0; i < headers.length; i++) {
			if (headers[i].equals("Stripped.Sequence")) {
				seqColIdx = i;
			} else if (headers[i].endsWith(".mzML")) {
				intensityColIndices.add(i);
			}
		}
		
		if (seqColIdx < 0) {
			reader.close();
			throw new IOException("Stripped.Sequence column not found");
		}
		if (intensityColIndices.isEmpty()) {
			reader.close();
			throw new IOException("No .mzML intensity columns found");
		}
		
		System.out.println("Found " + intensityColIndices.size() + " sample columns (.mzML)");
		
		// Store peptides with their intensities: peptide -> sample -> intensity
		Map<String, Map<String, Double>> peptideIntensities = new HashMap<>();
		int lineCount = 0;
		while ((headerLine = reader.readLine()) != null) {
			lineCount++;
			if (lineCount % 50000 == 0) {
				System.out.println("Read " + lineCount + " rows from matrix...");
			}
			
			String[] parts = headerLine.split("\t");
			if (parts.length <= seqColIdx) {
				continue;
			}
			
			String peptide = parts[seqColIdx].trim();
			if (peptide.isEmpty()) {
				continue;
			}
			
			// Collect intensities for this peptide
			Map<String, Double> intensities = new HashMap<>();
			for (int colIdx : intensityColIndices) {
				if (colIdx < parts.length) {
					try {
						double intensity = Double.parseDouble(parts[colIdx].trim());
						if (intensity > 0) {
							intensities.put(headers[colIdx], intensity);
						}
					} catch (NumberFormatException e) {
						// Skip invalid numbers
					}
				}
			}
			
			if (!intensities.isEmpty()) {
				peptideIntensities.put(peptide, intensities);
			}
		}
		reader.close();
		System.out.println("Loaded " + peptideIntensities.size() + " peptides with intensity data");
		
		System.out.println("Step 3: First pass - mapping peptides to genomes and identifying distinct peptides...");
		// Map peptides to genomes (only for peptides in experimental data)
		// peptide -> set of genomes (for ALL proteins)
		Map<String, Set<String>> peptideToGenomes = new ConcurrentHashMap<>();
		// peptide -> (genome -> category) for marker-4 proteins
		Map<String, Map<String, String>> peptideToGenomeCategory = new ConcurrentHashMap<>();
		
		ArrayList<File> peptideFiles = listPeptideFiles(predictedFolder);
		System.out.println("Processing " + peptideFiles.size() + " predicted TSV files with " + numThreads + " threads...");
		
		ExecutorService executor = Executors.newFixedThreadPool(numThreads);
		ArrayList<Future<?>> futures = new ArrayList<>();
		AtomicInteger processed = new AtomicInteger(0);
		
		// First pass: map peptides to genomes
		for (File pepFile : peptideFiles) {
			futures.add(executor.submit(() -> {
				String genomeId = genomeIdFromFile(pepFile.getName());
				if (genomeId == null) {
					return;
				}
				
				try {
					BufferedReader fileReader = new BufferedReader(new FileReader(pepFile));
					String line = fileReader.readLine(); // header
					while ((line = fileReader.readLine()) != null) {
						String[] parts = line.split("\t");
						if (parts.length < 2) {
							continue;
						}
						String proteinId = parts[0].trim();
						String peptide = parts[1].trim();
						
						// Only process peptides that exist in experimental data
						if (!peptideIntensities.containsKey(peptide)) {
							continue;
						}
						
						// Track ALL peptides for total abundance calculation
						peptideToGenomes.computeIfAbsent(peptide, k -> ConcurrentHashMap.newKeySet()).add(genomeId);
						
						// Track marker-4 peptides specifically
						String category = proteinToCategory.get(proteinId);
						if (category != null) {
							peptideToGenomeCategory.computeIfAbsent(peptide, k -> new ConcurrentHashMap<>()).put(genomeId, category);
						}
					}
					fileReader.close();
					
					int count = processed.incrementAndGet();
					if (count % 100 == 0) {
						System.out.println("Processed " + count + " / " + peptideFiles.size() + " files...");
						if (count % 500 == 0) {
							System.gc();
						}
					}
				} catch (IOException e) {
					System.err.println("Error reading " + pepFile.getName() + ": " + e.getMessage());
				}
			}));
		}
		
		// Wait for all tasks to complete
		for (Future<?> future : futures) {
			try {
				future.get();
			} catch (Exception e) {
				System.err.println("Error in parallel processing: " + e.getMessage());
			}
		}
		executor.shutdown();
		System.out.println("First pass completed. Found " + peptideToGenomeCategory.size() + " marker peptides.");
		
		// Identify distinct peptides (appear in only one genome)
		Set<String> distinctPeptides = new HashSet<>();
		for (Map.Entry<String, Set<String>> entry : peptideToGenomes.entrySet()) {
			if (entry.getValue() != null && entry.getValue().size() == 1) {
				distinctPeptides.add(entry.getKey());
			}
		}
		System.out.println("Found " + distinctPeptides.size() + " distinct peptides (unique to one genome)");
		
		System.out.println("Step 4: Calculating abundances (all peptides and distinct peptides)...");
		// Initialize abundance maps
		Map<String, Map<String, Map<String, Double>>> abundancesAll = new HashMap<>();
		Map<String, Map<String, Map<String, Double>>> abundancesDistinct = new HashMap<>();
		Map<String, Map<String, Double>> totalGenomeAbundances = new HashMap<>();
		
		for (String sampleName : headers) {
			if (!sampleName.endsWith(".mzML")) {
				continue;
			}
			abundancesAll.put(sampleName, new HashMap<>());
			abundancesDistinct.put(sampleName, new HashMap<>());
			totalGenomeAbundances.put(sampleName, new HashMap<>());
		}
		
		// Calculate abundances from peptide intensities
		for (Map.Entry<String, Map<String, Double>> peptideEntry : peptideIntensities.entrySet()) {
			String peptide = peptideEntry.getKey();
			Map<String, Double> intensities = peptideEntry.getValue();
			
			Set<String> genomes = peptideToGenomes.get(peptide);
			Map<String, String> genomeCategory = peptideToGenomeCategory.get(peptide);
			boolean isDistinct = distinctPeptides.contains(peptide);
			
			for (Map.Entry<String, Double> intensityEntry : intensities.entrySet()) {
				String sample = intensityEntry.getKey();
				double intensity = intensityEntry.getValue();
				
				// Add to total genome abundances (all proteins)
				if (genomes != null && !genomes.isEmpty()) {
					Map<String, Double> sampleTotals = totalGenomeAbundances.get(sample);
					if (sampleTotals != null) {
						for (String genomeId : genomes) {
							sampleTotals.put(genomeId, sampleTotals.getOrDefault(genomeId, 0.0) + intensity);
						}
					}
				}
				
				// Add to marker-4 category abundances (all peptides)
				if (genomeCategory != null && !genomeCategory.isEmpty()) {
					Map<String, Map<String, Double>> sampleAbundances = abundancesAll.get(sample);
					if (sampleAbundances != null) {
						for (Map.Entry<String, String> gcEntry : genomeCategory.entrySet()) {
							String genomeId = gcEntry.getKey();
							String category = gcEntry.getValue();
							
							sampleAbundances.putIfAbsent(genomeId, new HashMap<>());
							Map<String, Double> catAbundances = sampleAbundances.get(genomeId);
							catAbundances.put(category, catAbundances.getOrDefault(category, 0.0) + intensity);
						}
					}
				}
				
				// Add to marker-4 category abundances (distinct peptides only)
				if (isDistinct && genomeCategory != null && !genomeCategory.isEmpty()) {
					Set<String> distinctGenomes = peptideToGenomes.get(peptide);
					if (distinctGenomes != null && distinctGenomes.size() == 1) {
						String distinctGenomeId = distinctGenomes.iterator().next();
						String category = genomeCategory.get(distinctGenomeId);
						if (category != null) {
							Map<String, Map<String, Double>> sampleAbundances = abundancesDistinct.get(sample);
							if (sampleAbundances != null) {
								sampleAbundances.putIfAbsent(distinctGenomeId, new HashMap<>());
								Map<String, Double> catAbundances = sampleAbundances.get(distinctGenomeId);
								catAbundances.put(category, catAbundances.getOrDefault(category, 0.0) + intensity);
							}
						}
					}
				}
			}
		}
		
		System.out.println("Step 5: Calculating ratios and writing output...");
		// Write output with abundance ratios (both all peptides and distinct peptides)
		PrintWriter writer = new PrintWriter(outputFile);
		writer.println("Sample\tGenomeId\tCategory\tAbundance_All\tAbundance_Distinct\t" +
			"Ratio_All\tRatio_Distinct\tTotal_Genome_Abundance");
		
		String[] categories = {"Anchor", "HighRes", "Abundance", "Stability"};
		
		// Get all unique genomes across both all and distinct
		Set<String> allGenomes = new HashSet<>();
		for (Map<String, Map<String, Double>> sampleAbundances : abundancesAll.values()) {
			allGenomes.addAll(sampleAbundances.keySet());
		}
		for (Map<String, Map<String, Double>> sampleAbundances : abundancesDistinct.values()) {
			allGenomes.addAll(sampleAbundances.keySet());
		}
		
		for (String sample : abundancesAll.keySet()) {
			Map<String, Map<String, Double>> sampleAbundancesAll = abundancesAll.get(sample);
			Map<String, Map<String, Double>> sampleAbundancesDistinct = abundancesDistinct.get(sample);
			Map<String, Double> sampleTotals = totalGenomeAbundances.get(sample);
			
			// For each genome, calculate ratios
			for (String genomeId : allGenomes) {
				Map<String, Double> catAbundancesAll = sampleAbundancesAll.getOrDefault(genomeId, new HashMap<>());
				Map<String, Double> catAbundancesDistinct = sampleAbundancesDistinct.getOrDefault(genomeId, new HashMap<>());
				
				// Get total genome abundance (all proteins, not just marker-4)
				double totalGenomeAbundance = sampleTotals.getOrDefault(genomeId, 0.0);
				
				// Write all categories for this genome with ratios
				for (String category : categories) {
					double abundanceAll = catAbundancesAll.getOrDefault(category, 0.0);
					double abundanceDistinct = catAbundancesDistinct.getOrDefault(category, 0.0);
					
					// Calculate ratios: category abundance / total genome abundance (all proteins)
					double ratioAll = (totalGenomeAbundance > 0) ? (abundanceAll / totalGenomeAbundance) : 0.0;
					double ratioDistinct = (totalGenomeAbundance > 0) ? (abundanceDistinct / totalGenomeAbundance) : 0.0;
					
					writer.printf("%s\t%s\t%s\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f%n", 
						sample, genomeId, category, abundanceAll, abundanceDistinct, 
						ratioAll, ratioDistinct, totalGenomeAbundance);
				}
			}
		}
		
		writer.close();
		
		// Write summary: marker-4 percentage of total sample abundance
		PrintWriter summaryWriter = new PrintWriter(new File(outputFile.getParent(), 
			outputFile.getName().replace(".tsv", "_summary.tsv")));
		summaryWriter.println("Sample\tTotal_Marker4_Abundance_All\tTotal_Marker4_Abundance_Distinct\t" +
			"Total_All_Proteins_Abundance\tMarker4_Percentage_All\tMarker4_Percentage_Distinct");
		
		for (String sample : abundancesAll.keySet()) {
			// Calculate total marker-4 abundance for this sample (all peptides)
			double marker4TotalAll = 0;
			for (Map<String, Double> genomeAbundances : abundancesAll.get(sample).values()) {
				for (Double ab : genomeAbundances.values()) {
					marker4TotalAll += ab;
				}
			}
			
			// Calculate total marker-4 abundance for this sample (distinct peptides)
			double marker4TotalDistinct = 0;
			for (Map<String, Double> genomeAbundances : abundancesDistinct.get(sample).values()) {
				for (Double ab : genomeAbundances.values()) {
					marker4TotalDistinct += ab;
				}
			}
			
			// Calculate total all-proteins abundance for this sample
			double allProteinsTotal = totalGenomeAbundances.get(sample).values().stream()
				.mapToDouble(Double::doubleValue).sum();
			
			double percentageAll = (allProteinsTotal > 0) ? (100.0 * marker4TotalAll / allProteinsTotal) : 0.0;
			double percentageDistinct = (allProteinsTotal > 0) ? (100.0 * marker4TotalDistinct / allProteinsTotal) : 0.0;
			summaryWriter.printf("%s\t%.6f\t%.6f\t%.6f\t%.4f%%\t%.4f%%%n", 
				sample, marker4TotalAll, marker4TotalDistinct, allProteinsTotal, percentageAll, percentageDistinct);
		}
		summaryWriter.close();
		System.out.println("Summary file written: " + outputFile.getName().replace(".tsv", "_summary.tsv"));
		
		writer.close();
		System.out.println("Protein abundances written to: " + outputFile.getPath());
	}
	
	/**
	 * Convenience method using default thread count (CPU cores).
	 */
	public void calculateProteinAbundances(File matrixFile, File predictedFolder, File faaFolder, 
			File outputFile) throws IOException {
		calculateProteinAbundances(matrixFile, predictedFolder, faaFolder, outputFile, 0);
	}
	
	public static void main(String[] args) throws IOException {
		HapExtractor he = new HapExtractor();
		
		String dbFolder = "C:\\Users\\daz25hak\\MetaProAI\\database\\MAG\\human-gut\\v2.0.2\\original_db";
		String outputFolder = "C:\\Users\\daz25hak\\MetaPilot\\MetaPilot\\paper\\";
		
		// Original HAP method (ribosomal + elongation, strict)
//		he.getStatistics(new File(dbFolder), 
//				new File(outputFolder + "humangut_hap_strict.csv"));
		
		// 4-category method statistics (matching Python)
//		he.getCategoryStatistics(new File(dbFolder), 
//				new File(outputFolder + "humangut_4category.csv"));
		
		// Extract proteins from 4 categories to combined FASTA
		he.extractFourCategories(new File(dbFolder), 
				new File("C:\\Users\\daz25hak\\MetaProAI\\database\\MAG\\human-gut\\v2.0.2\\" + "4markers.fasta"));
		
		// Example: summarize peptide sharing for 4-category proteins
//		he.summarizePeptideSharing(new File(dbFolder),
//				new File(outputFolder),
//				new File(outputFolder + "humangut_4category_peptide_sharing.csv"));
		
//		he.generateProteinPeptideDistribution(
//				new File("C:\\Users\\daz25hak\\MetaProAI\\database\\MAG\\human-gut\\v2.0.2\\original_db"),
//				new File("C:\\Users\\daz25hak\\MetaProAI\\database\\MAG\\human-gut\\v2.0.2\\predicted"),
//				new File("C:\\Users\\daz25hak\\MetaPilot\\MetaPilot\\paper\\protein_peptide_distribution.tsv"));
		
//		he.generateDistinctPeptideDistribution(
//			    new File("C:\\Users\\daz25hak\\MetaProAI\\database\\MAG\\human-gut\\v2.0.2\\original_db"),
//			    new File("C:\\Users\\daz25hak\\MetaProAI\\database\\MAG\\human-gut\\v2.0.2\\predicted"),
//			    new File("C:\\Users\\daz25hak\\MetaPilot\\MetaPilot\\paper\\distinct_peptide_distribution.tsv")
//			);
		
//		he.calculateProteinAbundances(
//			    new File("C:\\Users\\daz25hak\\MetaPilot\\MetaPilot\\paper\\O480_60min_DIA_IGC_humanswiss.pr_matrix.tsv"),
//			    new File("C:\\Users\\daz25hak\\MetaProAI\\database\\MAG\\human-gut\\v2.0.2\\predicted"),
//			    new File("C:\\Users\\daz25hak\\MetaProAI\\database\\MAG\\human-gut\\v2.0.2\\original_db"),
//			    new File("C:\\Users\\daz25hak\\MetaPilot\\MetaPilot\\paper\\O480_60min_DIA_IGC_humanswiss_protein_abundances.tsv"),
//			    20  // Optional: number of threads (0 = auto-detect CPU cores)
//			);
	}
}

