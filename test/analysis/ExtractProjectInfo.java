import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * Extract complete information for selected projects from pride_datasets.json
 * Outputs to a single JSON file for further analysis
 */
public class ExtractProjectInfo {
	
	public static void main(String[] args) {
		// Try multiple possible paths
		File jsonFile = null;
		File projectsFile = null;
		
		String[] jsonPaths = {
			"pride_datasets.json",
			"./pride_datasets.json",
			"test/analysis/pride_datasets.json",
			"analysis/pride_datasets.json"
		};
		
		String[] projectPaths = {
			"projects.txt",
			"./projects.txt",
			"test/analysis/projects.txt",
			"analysis/projects.txt"
		};
		
		// Find JSON file
		for (String path : jsonPaths) {
			File f = new File(path);
			if (f.exists()) {
				jsonFile = f.getAbsoluteFile();
				break;
			}
		}
		
		// Find projects file
		for (String path : projectPaths) {
			File f = new File(path);
			if (f.exists()) {
				projectsFile = f.getAbsoluteFile();
				break;
			}
		}
		
		// If JSON found but projects not, try same directory as JSON
		if (jsonFile != null && jsonFile.exists() && (projectsFile == null || !projectsFile.exists())) {
			File parent = jsonFile.getParentFile();
			if (parent != null) {
				File testFile = new File(parent, "projects.txt");
				if (testFile.exists()) {
					projectsFile = testFile.getAbsoluteFile();
				}
			}
		}
		
		if (jsonFile == null || !jsonFile.exists()) {
			System.err.println("Error: pride_datasets.json not found");
			return;
		}
		
		if (projectsFile == null || !projectsFile.exists()) {
			System.err.println("Error: projects.txt not found");
			return;
		}
		
		File outputFile = new File(jsonFile.getParentFile(), "selected_projects_info.json");
		
		System.out.println("Reading JSON from: " + jsonFile.getAbsolutePath());
		System.out.println("Reading projects from: " + projectsFile.getAbsolutePath());
		System.out.println("Writing to: " + outputFile.getAbsolutePath());
		System.out.println("This may take a while for large JSON files...");
		
		try {
			// Read project list
			Set<String> projectSet = readProjects(projectsFile);
			System.out.println("Found " + projectSet.size() + " projects in projects.txt");
			
			// Extract matching project information
			List<String> projectJsonObjects = extractProjectInfo(jsonFile, projectSet);
			System.out.println("Extracted " + projectJsonObjects.size() + " matching projects");
			
			// Write to output file
			PrintWriter writer = new PrintWriter(outputFile);
			writer.println("[");
			for (int i = 0; i < projectJsonObjects.size(); i++) {
				writer.print(projectJsonObjects.get(i));
				if (i < projectJsonObjects.size() - 1) {
					writer.println(",");
				} else {
					writer.println();
				}
			}
			writer.println("]");
			writer.close();
			
			System.out.println("\nExtraction complete!");
			System.out.println("Output written to: " + outputFile.getAbsolutePath());
			System.out.println("You can now analyze this file with your preferred tools.");
			
		} catch (Exception e) {
			System.err.println("Error: " + e.getMessage());
			e.printStackTrace();
		}
	}
	
	private static Set<String> readProjects(File file) throws IOException {
		Set<String> projects = new HashSet<>();
		Pattern pxdPattern = Pattern.compile("(PXD\\d+)");
		
		BufferedReader reader = new BufferedReader(new FileReader(file));
		String line;
		while ((line = reader.readLine()) != null) {
			line = line.trim();
			if (line.isEmpty() || line.startsWith("#")) {
				continue;
			}
			
			Matcher matcher = pxdPattern.matcher(line);
			if (matcher.find()) {
				projects.add(matcher.group(1));
			}
		}
		reader.close();
		
		return projects;
	}
	
	/**
	 * Extract JSON objects for matching projects
	 */
	private static List<String> extractProjectInfo(File jsonFile, Set<String> projectSet) throws IOException {
		List<String> projectObjects = new ArrayList<>();
		
		BufferedReader reader = new BufferedReader(new FileReader(jsonFile));
		StringBuilder currentObject = new StringBuilder();
		int braceCount = 0;
		boolean inString = false;
		boolean escapeNext = false;
		
		String line;
		int lineNum = 0;
		int foundCount = 0;
		
		while ((line = reader.readLine()) != null) {
			lineNum++;
			if (lineNum % 10000 == 0) {
				System.out.println("Processing line " + lineNum + "... (found " + foundCount + " projects)");
			}
			
			for (char c : line.toCharArray()) {
				if (escapeNext) {
					escapeNext = false;
					currentObject.append(c);
					continue;
				}
				
				if (c == '\\') {
					escapeNext = true;
					currentObject.append(c);
					continue;
				}
				
				if (c == '"') {
					inString = !inString;
					currentObject.append(c);
					continue;
				}
				
				if (!inString) {
					if (c == '{') {
						braceCount++;
						currentObject.append(c);
					} else if (c == '}') {
						braceCount--;
						currentObject.append(c);
						
						if (braceCount == 0) {
							// Complete object found
							String jsonObject = currentObject.toString();
							
							// Check if this project is in our list
							Pattern accessionPattern = Pattern.compile("\"accession\"\\s*:\\s*\"(PXD\\d+)\"");
							Matcher matcher = accessionPattern.matcher(jsonObject);
							if (matcher.find()) {
								String accession = matcher.group(1);
								if (projectSet.contains(accession)) {
									projectObjects.add(jsonObject);
									foundCount++;
									if (foundCount % 10 == 0) {
										System.out.println("  Found project " + foundCount + ": " + accession);
									}
								}
							}
							
							currentObject.setLength(0);
						}
					} else {
						currentObject.append(c);
					}
				} else {
					currentObject.append(c);
				}
			}
			
			currentObject.append('\n');
		}
		
		reader.close();
		return projectObjects;
	}
}

