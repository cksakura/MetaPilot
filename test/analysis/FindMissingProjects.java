import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.HashSet;
import java.util.Set;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.TreeSet;

/**
 * Find projects from projects.txt that don't have any .raw files in raw_files_list.txt
 */
public class FindMissingProjects {
	
	public static void main(String[] args) {
		// Try multiple possible paths - use absolute paths to avoid duplication
		File projectsFile = null;
		File rawFilesFile = null;
		File outputFile = null;
		
		String[] possiblePaths = {
			"projects.txt",
			"./projects.txt",
			"test/analysis/projects.txt",
			"analysis/projects.txt"
		};
		
		for (String path : possiblePaths) {
			File f = new File(path);
			if (f.exists()) {
				projectsFile = f.getAbsoluteFile();
				// Use the same directory as projects.txt for other files
				File baseDir = projectsFile.getParentFile();
				if (baseDir != null) {
					rawFilesFile = new File(baseDir, "raw_files_list.txt").getAbsoluteFile();
					outputFile = new File(baseDir, "missing_projects.txt").getAbsoluteFile();
				} else {
					// If no parent, use current directory
					rawFilesFile = new File("raw_files_list.txt").getAbsoluteFile();
					outputFile = new File("missing_projects.txt").getAbsoluteFile();
				}
				break;
			}
		}
		
		if (projectsFile == null || !projectsFile.exists()) {
			System.err.println("Error: projects.txt not found");
			System.err.println("Searched in current directory and common paths");
			return;
		}
		
		if (rawFilesFile == null || !rawFilesFile.exists()) {
			System.err.println("Error: raw_files_list.txt not found");
			System.err.println("Expected at: " + rawFilesFile.getAbsolutePath());
			return;
		}
		
		System.out.println("Reading projects from: " + projectsFile.getAbsolutePath());
		System.out.println("Reading raw files from: " + rawFilesFile.getAbsolutePath());
		System.out.println("Writing missing projects to: " + outputFile.getAbsolutePath());
		
		try {
			// Pattern to extract PXD ID from URL
			Pattern pxdPattern = Pattern.compile("(PXD\\d+)");
			
			// Read all projects from projects.txt
			Set<String> allProjects = new TreeSet<>();
			Set<String> projectLinks = new HashSet<>();
			
			BufferedReader reader = new BufferedReader(new FileReader(projectsFile));
			String line;
			while ((line = reader.readLine()) != null) {
				line = line.trim();
				if (line.isEmpty() || line.startsWith("#")) {
					continue;
				}
				
				Matcher matcher = pxdPattern.matcher(line);
				if (matcher.find()) {
					String pxdId = matcher.group(1);
					allProjects.add(pxdId);
					projectLinks.add(line);
				}
			}
			reader.close();
			
			System.out.println("\nTotal projects in projects.txt: " + allProjects.size());
			
			// Read projects that have .raw files from raw_files_list.txt
			Set<String> projectsWithRawFiles = new HashSet<>();
			
			reader = new BufferedReader(new FileReader(rawFilesFile));
			line = reader.readLine(); // Skip header
			while ((line = reader.readLine()) != null) {
				line = line.trim();
				if (line.isEmpty()) {
					continue;
				}
				
				String[] parts = line.split("\t");
				if (parts.length < 1) {
					continue;
				}
				
				String ftpLink = parts[0];
				Matcher matcher = pxdPattern.matcher(ftpLink);
				if (matcher.find()) {
					String pxdId = matcher.group(1);
					projectsWithRawFiles.add(pxdId);
				}
			}
			reader.close();
			
			System.out.println("Projects with .raw files: " + projectsWithRawFiles.size());
			
			// Find missing projects
			Set<String> missingProjects = new TreeSet<>();
			for (String pxdId : allProjects) {
				if (!projectsWithRawFiles.contains(pxdId)) {
					missingProjects.add(pxdId);
				}
			}
			
			System.out.println("Missing projects (no .raw files): " + missingProjects.size());
			
			// Write missing projects with their full URLs
			PrintWriter writer = new PrintWriter(outputFile);
			writer.println("Missing Projects (No .raw files found)");
			writer.println("=".repeat(80));
			writer.println();
			writer.println("PXD_ID\tFTP_Link");
			writer.println("-".repeat(80));
			
			for (String pxdId : missingProjects) {
				// Find the corresponding URL from projectLinks
				String url = null;
				for (String link : projectLinks) {
					if (link.contains(pxdId)) {
						url = link;
						break;
					}
				}
				
				if (url != null) {
					writer.println(pxdId + "\t" + url);
				} else {
					writer.println(pxdId + "\t(URL not found)");
				}
			}
			
			writer.println();
			writer.println("=".repeat(80));
			writer.println("Total missing: " + missingProjects.size());
			writer.close();
			
			// Print summary to console
			System.out.println("\nMissing Projects:");
			System.out.println("-".repeat(80));
			for (String pxdId : missingProjects) {
				// Find the corresponding URL
				String url = null;
				for (String link : projectLinks) {
					if (link.contains(pxdId)) {
						url = link;
						break;
					}
				}
				if (url != null) {
					System.out.println(pxdId + " - " + url);
				} else {
					System.out.println(pxdId + " - (URL not found)");
				}
			}
			
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

