import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.TreeMap;

/**
 * Analyze raw_files_list.txt to generate statistics per year:
 * - Number of projects per year
 * - Number of raw files per year
 */
public class ProjectStatistics {
	
	public static void main(String[] args) {
		File inputFile = new File("test/analysis/raw_files_list.txt");
		File outputFile = new File("test/analysis/project_statistics.txt");
		
		// Try multiple possible paths
		if (!inputFile.exists()) {
			String[] possiblePaths = {
				"raw_files_list.txt",
				"./raw_files_list.txt",
				"analysis/raw_files_list.txt"
			};
			for (String path : possiblePaths) {
				File f = new File(path);
				if (f.exists()) {
					inputFile = f;
					outputFile = new File(f.getParent(), "project_statistics.txt");
					break;
				}
			}
		}
		
		if (!inputFile.exists()) {
			System.err.println("Error: raw_files_list.txt not found");
			System.err.println("Tried: " + inputFile.getAbsolutePath());
			return;
		}
		
		System.out.println("Reading from: " + inputFile.getAbsolutePath());
		System.out.println("Writing to: " + outputFile.getAbsolutePath());
		
		try {
			// Parse the file and collect statistics
			Map<Integer, Set<String>> projectsPerYear = new TreeMap<>(); // Year -> Set of PXD IDs
			Map<Integer, Integer> rawFilesPerYear = new HashMap<>(); // Year -> Count of raw files
			
			// Pattern to extract year from URL: /archive/YYYY/
			Pattern yearPattern = Pattern.compile("/archive/(\\d{4})/");
			// Pattern to extract PXD ID: PXD followed by digits
			Pattern pxdPattern = Pattern.compile("(PXD\\d+)");
			
			BufferedReader reader = new BufferedReader(new FileReader(inputFile));
			String line = reader.readLine(); // Skip header
			
			int totalLines = 0;
			while ((line = reader.readLine()) != null) {
				totalLines++;
				line = line.trim();
				if (line.isEmpty()) {
					continue;
				}
				
				String[] parts = line.split("\t");
				if (parts.length < 2) {
					continue;
				}
				
				String ftpLink = parts[0];
				String rawFileName = parts[1];
				
				// Extract year from FTP link
				Matcher yearMatcher = yearPattern.matcher(ftpLink);
				if (yearMatcher.find()) {
					int year = Integer.parseInt(yearMatcher.group(1));
					
					// Extract PXD ID
					Matcher pxdMatcher = pxdPattern.matcher(ftpLink);
					if (pxdMatcher.find()) {
						String pxdId = pxdMatcher.group(1);
						
						// Add project to this year
						projectsPerYear.putIfAbsent(year, new HashSet<>());
						projectsPerYear.get(year).add(pxdId);
						
						// Count raw file for this year
						rawFilesPerYear.put(year, rawFilesPerYear.getOrDefault(year, 0) + 1);
					}
				}
			}
			reader.close();
			
			// Write statistics
			PrintWriter writer = new PrintWriter(outputFile);
			writer.println("Year\tProjects\tRaw_Files\tAvg_Raw_Files_per_Project");
			writer.println(repeatString("=", 60));
			
			int totalProjects = 0;
			int totalRawFiles = 0;
			
			for (Map.Entry<Integer, Set<String>> entry : projectsPerYear.entrySet()) {
				int year = entry.getKey();
				int projectCount = entry.getValue().size();
				int rawFileCount = rawFilesPerYear.getOrDefault(year, 0);
				double avgRawFiles = projectCount > 0 ? (double) rawFileCount / projectCount : 0.0;
				
				writer.printf("%d\t%d\t%d\t%.2f%n", year, projectCount, rawFileCount, avgRawFiles);
				
				totalProjects += projectCount;
				totalRawFiles += rawFileCount;
			}
			
			writer.println(repeatString("=", 60));
			writer.printf("TOTAL\t%d\t%d\t%.2f%n", totalProjects, totalRawFiles, 
				totalProjects > 0 ? (double) totalRawFiles / totalProjects : 0.0);
			writer.close();
			
			// Print summary to console
			System.out.println("\nStatistics Summary:");
			System.out.println(repeatString("=", 60));
			System.out.printf("%-10s %-15s %-15s %-25s%n", "Year", "Projects", "Raw Files", "Avg Raw/Project");
			System.out.println(repeatString("=", 60));
			
			for (Map.Entry<Integer, Set<String>> entry : projectsPerYear.entrySet()) {
				int year = entry.getKey();
				int projectCount = entry.getValue().size();
				int rawFileCount = rawFilesPerYear.getOrDefault(year, 0);
				double avgRawFiles = projectCount > 0 ? (double) rawFileCount / projectCount : 0.0;
				
				System.out.printf("%-10d %-15d %-15d %-25.2f%n", year, projectCount, rawFileCount, avgRawFiles);
			}
			
			System.out.println(repeatString("=", 60));
			System.out.printf("%-10s %-15d %-15d %-25.2f%n", "TOTAL", totalProjects, totalRawFiles,
				totalProjects > 0 ? (double) totalRawFiles / totalProjects : 0.0);
			System.out.println("\nTotal lines processed: " + totalLines);
			System.out.println("\nResults written to: " + outputFile.getAbsolutePath());
			
		} catch (Exception e) {
			System.err.println("Error: " + e.getMessage());
			e.printStackTrace();
		}
	}
	
	private static String repeatString(String str, int count) {
		StringBuilder sb = new StringBuilder();
		for (int i = 0; i < count; i++) {
			sb.append(str);
		}
		return sb.toString();
	}
}

