import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;

/**
 * Combine multiple result files from Figure1D, Figure1E, and Figure1F folders
 * into single combined files for visualization
 */
public class CombineResults {
	
	public static void main(String[] args) {
		// Get the directory where the script is located
		File currentDir = new File(System.getProperty("user.dir"));
		File baseDir;
		
		// Check if we're in the test/analysis directory or need to navigate
		if (new File(currentDir, "Figure1D").exists()) {
			baseDir = currentDir;
		} else if (new File(currentDir, "test/analysis/Figure1D").exists()) {
			baseDir = new File(currentDir, "test/analysis");
		} else {
			// Try relative to script location
			baseDir = new File(".");
			if (!new File(baseDir, "Figure1D").exists()) {
				System.err.println("Error: Cannot find Figure1D folder. Please run from test/analysis directory.");
				return;
			}
		}
		
		// Combine files in each folder
		combineFolder(baseDir, "Figure1D", "combined_marker_abundances.tsv");
		combineFolder(baseDir, "Figure1E", "combined_marker_abundances_unique.tsv");
		combineFolder(baseDir, "Figure1F", "combined_peptide_detectability.tsv");
		
		System.out.println("\n=== Combination Complete ===");
		System.out.println("Output files:");
		System.out.println("  - " + new File(baseDir, "combined_marker_abundances.tsv").getAbsolutePath());
		System.out.println("  - " + new File(baseDir, "combined_marker_abundances_unique.tsv").getAbsolutePath());
		System.out.println("  - " + new File(baseDir, "combined_peptide_detectability.tsv").getAbsolutePath());
	}
	
	/**
	 * Combine all files in a folder into a single file
	 */
	private static void combineFolder(File baseDir, String folderName, String outputFileName) {
		File folder = new File(baseDir, folderName);
		
		if (!folder.exists() || !folder.isDirectory()) {
			System.err.println("Warning: Folder does not exist: " + folder.getPath());
			return;
		}
		
		File[] files = folder.listFiles();
		if (files == null || files.length == 0) {
			System.err.println("Warning: No files found in: " + folder.getPath());
			return;
		}
		
		// Filter to only TSV files and sort by name
		ArrayList<File> tsvFiles = new ArrayList<>();
		for (File f : files) {
			if (f.isFile() && (f.getName().endsWith(".tsv") || f.getName().endsWith(".txt"))) {
				tsvFiles.add(f);
			}
		}
		
		if (tsvFiles.isEmpty()) {
			System.err.println("Warning: No TSV/TXT files found in: " + folder.getPath());
			return;
		}
		
		// Sort files by name for consistent ordering
		tsvFiles.sort(Comparator.comparing(File::getName));
		
		System.out.println("Combining " + tsvFiles.size() + " files from " + folderName + "...");
		
		File outputFile = new File(baseDir, outputFileName);
		
		try {
			PrintWriter writer = new PrintWriter(outputFile);
			boolean isFirstFile = true;
			int totalLines = 0;
			
			for (File inputFile : tsvFiles) {
				System.out.println("  Processing: " + inputFile.getName());
				
				try (BufferedReader reader = new BufferedReader(new FileReader(inputFile))) {
					String headerLine = reader.readLine();
					
					if (headerLine == null) {
						System.err.println("    Warning: Empty file: " + inputFile.getName());
						continue;
					}
					
					// Write header only for the first file
					if (isFirstFile) {
						writer.println(headerLine);
						isFirstFile = false;
						totalLines++;
					}
					
					// Write all data lines
					String line;
					int lineCount = 0;
					while ((line = reader.readLine()) != null) {
						if (!line.trim().isEmpty()) {  // Skip empty lines
							writer.println(line);
							lineCount++;
							totalLines++;
						}
					}
					
					System.out.println("    Added " + lineCount + " data lines");
				} catch (IOException e) {
					System.err.println("    Error reading file " + inputFile.getName() + ": " + e.getMessage());
				}
			}
			
			writer.close();
			System.out.println("  Combined file created: " + outputFile.getName() + " (" + totalLines + " total lines)");
			
		} catch (IOException e) {
			System.err.println("Error writing combined file: " + e.getMessage());
			e.printStackTrace();
		}
	}
}

