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
 * Extract project information to CSV format for easy analysis in Excel/Python
 */
public class ExtractProjectInfoCSV {
	
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
		
		File outputFile = new File(jsonFile.getParentFile(), "selected_projects_info.csv");
		
		System.out.println("Reading JSON from: " + jsonFile.getAbsolutePath());
		System.out.println("Reading projects from: " + projectsFile.getAbsolutePath());
		System.out.println("Writing to: " + outputFile.getAbsolutePath());
		System.out.println("This may take a while for large JSON files...");
		
		try {
			// Read project list
			Set<String> projectSet = readProjects(projectsFile);
			System.out.println("Found " + projectSet.size() + " projects in projects.txt");
			
			// Extract matching project information
			List<ProjectInfo> projects = extractProjectInfo(jsonFile, projectSet);
			System.out.println("Extracted " + projects.size() + " matching projects");
			
			// Write to CSV
			PrintWriter writer = new PrintWriter(outputFile);
			// Write header
			writer.println("Accession,Title,Description,Keywords,Submission_Date,Publication_Date,DOI,License,Submission_Type,Submitter_Name,Submitter_Affiliation,Submitter_Email,PI_Name,PI_Affiliation,PI_Email,Instruments,Softwares");
			
			// Write data
			for (ProjectInfo info : projects) {
				writer.println(
					escapeCsv(info.accession) + "," +
					escapeCsv(info.title) + "," +
					escapeCsv(info.description) + "," +
					escapeCsv(String.join("; ", info.keywords)) + "," +
					escapeCsv(info.submissionDate) + "," +
					escapeCsv(info.publicationDate) + "," +
					escapeCsv(info.doi) + "," +
					escapeCsv(info.license) + "," +
					escapeCsv(info.submissionType) + "," +
					escapeCsv(info.submitterName) + "," +
					escapeCsv(info.submitterAffiliation) + "," +
					escapeCsv(info.submitterEmail) + "," +
					escapeCsv(info.piName) + "," +
					escapeCsv(info.piAffiliation) + "," +
					escapeCsv(info.piEmail) + "," +
					escapeCsv(String.join("; ", info.instruments)) + "," +
					escapeCsv(String.join("; ", info.softwares))
				);
			}
			
			writer.close();
			
			System.out.println("\nExtraction complete!");
			System.out.println("Output written to: " + outputFile.getAbsolutePath());
			System.out.println("You can now open this CSV file in Excel or analyze with Python/R.");
			
		} catch (Exception e) {
			System.err.println("Error: " + e.getMessage());
			e.printStackTrace();
		}
	}
	
	private static String escapeCsv(String value) {
		if (value == null || value.isEmpty()) {
			return "";
		}
		// Escape quotes and wrap in quotes if contains comma, quote, or newline
		value = value.replace("\"", "\"\"");
		if (value.contains(",") || value.contains("\"") || value.contains("\n") || value.contains("\r")) {
			return "\"" + value + "\"";
		}
		return value;
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
	 * Extract project information from JSON
	 */
	private static List<ProjectInfo> extractProjectInfo(File jsonFile, Set<String> projectSet) throws IOException {
		List<ProjectInfo> projects = new ArrayList<>();
		
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
							ProjectInfo info = parseProjectInfo(jsonObject, projectSet);
							if (info != null) {
								projects.add(info);
								foundCount++;
								if (foundCount % 10 == 0) {
									System.out.println("  Found project " + foundCount + ": " + info.accession);
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
		return projects;
	}
	
	/**
	 * Parse project information from JSON object
	 */
	private static ProjectInfo parseProjectInfo(String json, Set<String> projectSet) {
		ProjectInfo info = new ProjectInfo();
		
		// Extract accession
		Pattern accessionPattern = Pattern.compile("\"accession\"\\s*:\\s*\"(PXD\\d+)\"");
		Matcher matcher = accessionPattern.matcher(json);
		if (!matcher.find()) {
			return null;
		}
		info.accession = matcher.group(1);
		
		// Only process if in our project list
		if (!projectSet.contains(info.accession)) {
			return null;
		}
		
		// Extract fields using patterns
		extractField(json, "\"title\"\\s*:\\s*\"([^\"]+)\"", (val) -> info.title = val);
		extractField(json, "\"projectDescription\"\\s*:\\s*\"([^\"]+)\"", (val) -> info.description = val);
		extractField(json, "\"submissionDate\"\\s*:\\s*\"([^\"]+)\"", (val) -> info.submissionDate = val);
		extractField(json, "\"publicationDate\"\\s*:\\s*\"([^\"]+)\"", (val) -> info.publicationDate = val);
		extractField(json, "\"doi\"\\s*:\\s*\"([^\"]+)\"", (val) -> info.doi = val);
		extractField(json, "\"license\"\\s*:\\s*\"([^\"]+)\"", (val) -> info.license = val);
		extractField(json, "\"submissionType\"\\s*:\\s*\"([^\"]+)\"", (val) -> info.submissionType = val);
		
		// Extract keywords
		Pattern keywordsPattern = Pattern.compile("\"keywords\"\\s*:\\s*\\[([^\\]]+)\\]");
		matcher = keywordsPattern.matcher(json);
		if (matcher.find()) {
			String keywordsStr = matcher.group(1);
			Pattern keywordPattern = Pattern.compile("\"([^\"]+)\"");
			Matcher keywordMatcher = keywordPattern.matcher(keywordsStr);
			while (keywordMatcher.find()) {
				info.keywords.add(keywordMatcher.group(1));
			}
		}
		
		// Extract submitter info (first one)
		Pattern submitterPattern = Pattern.compile("\"submitters\"\\s*:\\s*\\[\\s*\\{([^}]+)\\}");
		matcher = submitterPattern.matcher(json);
		if (matcher.find()) {
			String submitterJson = matcher.group(1);
			extractField(submitterJson, "\"name\"\\s*:\\s*\"([^\"]+)\"", (val) -> info.submitterName = val);
			extractField(submitterJson, "\"affiliation\"\\s*:\\s*\"([^\"]+)\"", (val) -> info.submitterAffiliation = val);
			extractField(submitterJson, "\"email\"\\s*:\\s*\"([^\"]+)\"", (val) -> info.submitterEmail = val);
		}
		
		// Extract PI info (first one)
		Pattern piPattern = Pattern.compile("\"labPIs\"\\s*:\\s*\\[\\s*\\{([^}]+)\\}");
		matcher = piPattern.matcher(json);
		if (matcher.find()) {
			String piJson = matcher.group(1);
			extractField(piJson, "\"name\"\\s*:\\s*\"([^\"]+)\"", (val) -> info.piName = val);
			extractField(piJson, "\"affiliation\"\\s*:\\s*\"([^\"]+)\"", (val) -> info.piAffiliation = val);
			extractField(piJson, "\"email\"\\s*:\\s*\"([^\"]+)\"", (val) -> info.piEmail = val);
		}
		
		// Extract instruments
		Pattern instrumentPattern = Pattern.compile("\"instruments\"\\s*:\\s*\\[([^\\]]+)\\]");
		matcher = instrumentPattern.matcher(json);
		if (matcher.find()) {
			String instrumentsStr = matcher.group(1);
			Pattern namePattern = Pattern.compile("\"name\"\\s*:\\s*\"([^\"]+)\"");
			Matcher nameMatcher = namePattern.matcher(instrumentsStr);
			while (nameMatcher.find()) {
				info.instruments.add(nameMatcher.group(1));
			}
		}
		
		// Extract softwares
		Pattern softwarePattern = Pattern.compile("\"softwares\"\\s*:\\s*\\[([^\\]]+)\\]");
		matcher = softwarePattern.matcher(json);
		if (matcher.find()) {
			String softwaresStr = matcher.group(1);
			Pattern namePattern = Pattern.compile("\"name\"\\s*:\\s*\"([^\"]+)\"");
			Matcher nameMatcher = namePattern.matcher(softwaresStr);
			while (nameMatcher.find()) {
				info.softwares.add(nameMatcher.group(1));
			}
		}
		
		return info;
	}
	
	private static void extractField(String json, String pattern, java.util.function.Consumer<String> setter) {
		Pattern p = Pattern.compile(pattern);
		Matcher m = p.matcher(json);
		if (m.find()) {
			setter.accept(m.group(1));
		}
	}
	
	/**
	 * Project information class
	 */
	private static class ProjectInfo {
		String accession = "";
		String title = "";
		String description = "";
		List<String> keywords = new ArrayList<>();
		String submissionDate = "";
		String publicationDate = "";
		String doi = "";
		String license = "";
		String submissionType = "";
		String submitterName = "";
		String submitterAffiliation = "";
		String submitterEmail = "";
		String piName = "";
		String piAffiliation = "";
		String piEmail = "";
		List<String> instruments = new ArrayList<>();
		List<String> softwares = new ArrayList<>();
	}
}

