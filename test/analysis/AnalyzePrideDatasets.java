import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * Analyze PRIDE datasets JSON file to:
 * 1. Match with projects.txt
 * 2. Extract keywords and create word cloud data
 * 3. Count disease-related projects
 * 4. Generate summary statistics
 */
public class AnalyzePrideDatasets {
	
	// Common disease-related keywords
	private static final String[] DISEASE_KEYWORDS = {
		"cancer", "tumor", "carcinoma", "malignancy", "oncology",
		"diabetes", "diabetic",
		"alzheimer", "dementia", "parkinson",
		"cardiovascular", "heart", "cardiac", "myocardial",
		"obesity", "obese",
		"inflammation", "inflammatory",
		"infection", "infectious", "bacterial", "viral", "pathogen",
		"disease", "disorder", "syndrome", "pathology",
		"covid", "sars", "coronavirus",
		"asthma", "copd",
		"arthritis", "rheumatoid",
		"liver", "hepatic", "cirrhosis",
		"kidney", "renal", "nephritis",
		"autoimmune", "immune",
		"depression", "anxiety", "mental",
		"stroke", "ischemia"
	};
	
	public static void main(String[] args) {
		// Try multiple possible paths - use absolute paths to avoid duplication
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
			System.err.println("Searched in:");
			for (String path : jsonPaths) {
				System.err.println("  " + new File(path).getAbsolutePath());
			}
			return;
		}
		
		if (projectsFile == null || !projectsFile.exists()) {
			System.err.println("Error: projects.txt not found");
			System.err.println("Searched in:");
			for (String path : projectPaths) {
				System.err.println("  " + new File(path).getAbsolutePath());
			}
			if (jsonFile.getParentFile() != null) {
				System.err.println("  " + new File(jsonFile.getParentFile(), "projects.txt").getAbsolutePath());
			}
			return;
		}
		
		System.out.println("Reading JSON from: " + jsonFile.getAbsolutePath());
		System.out.println("Reading projects from: " + projectsFile.getAbsolutePath());
		System.out.println("This may take a while for large files...");
		
		try {
			// Read project list
			Set<String> projectSet = readProjects(projectsFile);
			System.out.println("Found " + projectSet.size() + " projects in projects.txt");
			
			// Parse JSON and extract information
			Map<String, DatasetInfo> datasets = parseJsonFile(jsonFile, projectSet);
			System.out.println("Found " + datasets.size() + " matching datasets");
			
			// Generate statistics
			generateStatistics(datasets, jsonFile.getParentFile());
			
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
	 * Parse JSON file line by line (streaming approach for large files)
	 */
	private static Map<String, DatasetInfo> parseJsonFile(File file, Set<String> projectSet) throws IOException {
		Map<String, DatasetInfo> datasets = new HashMap<>();
		
		BufferedReader reader = new BufferedReader(new FileReader(file));
		StringBuilder currentObject = new StringBuilder();
		int braceCount = 0;
		boolean inString = false;
		boolean escapeNext = false;
		
		String line;
		int lineNum = 0;
		
		while ((line = reader.readLine()) != null) {
			lineNum++;
			if (lineNum % 10000 == 0) {
				System.out.println("Processing line " + lineNum + "...");
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
							DatasetInfo info = parseDatasetObject(jsonObject, projectSet);
							if (info != null) {
								datasets.put(info.accession, info);
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
		return datasets;
	}
	
	/**
	 * Parse a single dataset JSON object
	 */
	private static DatasetInfo parseDatasetObject(String json, Set<String> projectSet) {
		DatasetInfo info = new DatasetInfo();
		
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
		
		// Extract title
		Pattern titlePattern = Pattern.compile("\"title\"\\s*:\\s*\"([^\"]+)\"");
		matcher = titlePattern.matcher(json);
		if (matcher.find()) {
			info.title = matcher.group(1);
		}
		
		// Extract project description
		Pattern descPattern = Pattern.compile("\"projectDescription\"\\s*:\\s*\"([^\"]+)\"");
		matcher = descPattern.matcher(json);
		if (matcher.find()) {
			info.description = matcher.group(1);
		}
		
		// Extract keywords
		Pattern keywordsPattern = Pattern.compile("\"keywords\"\\s*:\\s*\\[([^\\]]+)\\]");
		matcher = keywordsPattern.matcher(json);
		if (matcher.find()) {
			String keywordsStr = matcher.group(1);
			// Parse individual keywords
			Pattern keywordPattern = Pattern.compile("\"([^\"]+)\"");
			Matcher keywordMatcher = keywordPattern.matcher(keywordsStr);
			while (keywordMatcher.find()) {
				info.keywords.add(keywordMatcher.group(1));
			}
		}
		
		// Check if disease-related
		info.isDiseaseRelated = checkDiseaseRelated(info);
		
		return info;
	}
	
	/**
	 * Check if dataset is disease-related
	 */
	private static boolean checkDiseaseRelated(DatasetInfo info) {
		String text = (info.title + " " + info.description + " " + String.join(" ", info.keywords)).toLowerCase();
		
		for (String keyword : DISEASE_KEYWORDS) {
			if (text.contains(keyword.toLowerCase())) {
				return true;
			}
		}
		
		return false;
	}
	
	/**
	 * Generate statistics and output files
	 */
	private static void generateStatistics(Map<String, DatasetInfo> datasets, File outputDir) throws IOException {
		// Count disease-related projects
		int diseaseCount = 0;
		for (DatasetInfo info : datasets.values()) {
			if (info.isDiseaseRelated) {
				diseaseCount++;
			}
		}
		
		// Collect all words for word cloud
		Map<String, Integer> wordCounts = new HashMap<>();
		for (DatasetInfo info : datasets.values()) {
			// Extract words from title, description, and keywords
			String text = info.title + " " + info.description + " " + String.join(" ", info.keywords);
			extractWords(text, wordCounts);
		}
		
		// Write word cloud data
		File wordCloudFile = new File(outputDir, "word_cloud_data.csv");
		PrintWriter writer = new PrintWriter(wordCloudFile);
		writer.println("word,count");
		
		// Sort by count and write top words
		final PrintWriter finalWriter = writer; // Make final for lambda
		wordCounts.entrySet().stream()
			.sorted((a, b) -> b.getValue().compareTo(a.getValue()))
			.limit(500) // Top 500 words
			.forEach(entry -> finalWriter.println(entry.getKey() + "," + entry.getValue()));
		
		writer.close();
		
		// Write summary statistics
		File summaryFile = new File(outputDir, "dataset_summary.txt");
		writer = new PrintWriter(summaryFile);
		writer.println("PRIDE Dataset Analysis Summary");
		writer.println("=".repeat(60));
		writer.println();
		writer.println("Total projects analyzed: " + datasets.size());
		writer.println("Disease-related projects: " + diseaseCount);
		writer.println("Non-disease projects: " + (datasets.size() - diseaseCount));
		writer.println("Percentage disease-related: " + String.format("%.2f%%", (double)diseaseCount / datasets.size() * 100));
		writer.println();
		writer.println("Word cloud data written to: word_cloud_data.csv");
		writer.println("(Top 500 words by frequency)");
		writer.close();
		
		// Write disease-related projects list
		File diseaseFile = new File(outputDir, "disease_related_projects.txt");
		writer = new PrintWriter(diseaseFile);
		writer.println("Disease-Related Projects");
		writer.println("=".repeat(60));
		writer.println();
		
		for (DatasetInfo info : datasets.values()) {
			if (info.isDiseaseRelated) {
				writer.println(info.accession + "\t" + info.title);
			}
		}
		writer.close();
		
		// Print summary to console
		System.out.println("\n" + "=".repeat(60));
		System.out.println("Analysis Complete!");
		System.out.println("=".repeat(60));
		System.out.println("Total projects: " + datasets.size());
		System.out.println("Disease-related: " + diseaseCount + " (" + String.format("%.2f%%", (double)diseaseCount / datasets.size() * 100) + ")");
		System.out.println("Non-disease: " + (datasets.size() - diseaseCount));
		System.out.println();
		System.out.println("Output files:");
		System.out.println("  - word_cloud_data.csv (for word cloud visualization)");
		System.out.println("  - dataset_summary.txt (summary statistics)");
		System.out.println("  - disease_related_projects.txt (list of disease projects)");
	}
	
	/**
	 * Extract words from text (simple tokenization)
	 */
	private static void extractWords(String text, Map<String, Integer> wordCounts) {
		// Remove special characters and split
		text = text.toLowerCase().replaceAll("[^a-z0-9\\s]", " ");
		String[] words = text.split("\\s+");
		
		// Common stop words to exclude
		Set<String> stopWords = new HashSet<>();
		stopWords.add("the"); stopWords.add("a"); stopWords.add("an"); stopWords.add("and"); stopWords.add("or");
		stopWords.add("but"); stopWords.add("in"); stopWords.add("on"); stopWords.add("at"); stopWords.add("to");
		stopWords.add("for"); stopWords.add("of"); stopWords.add("with"); stopWords.add("by"); stopWords.add("from");
		stopWords.add("is"); stopWords.add("are"); stopWords.add("was"); stopWords.add("were"); stopWords.add("be");
		stopWords.add("been"); stopWords.add("being"); stopWords.add("have"); stopWords.add("has"); stopWords.add("had");
		stopWords.add("do"); stopWords.add("does"); stopWords.add("did"); stopWords.add("will"); stopWords.add("would");
		stopWords.add("this"); stopWords.add("that"); stopWords.add("these"); stopWords.add("those");
		
		for (String word : words) {
			word = word.trim();
			if (word.length() > 2 && !stopWords.contains(word) && !word.matches("^\\d+$")) {
				wordCounts.put(word, wordCounts.getOrDefault(word, 0) + 1);
			}
		}
	}
	
	/**
	 * Dataset information class
	 */
	private static class DatasetInfo {
		String accession;
		String title = "";
		String description = "";
		List<String> keywords = new ArrayList<>();
		boolean isDiseaseRelated = false;
	}
	
	private static String repeatString(String str, int count) {
		StringBuilder sb = new StringBuilder();
		for (int i = 0; i < count; i++) {
			sb.append(str);
		}
		return sb.toString();
	}
}

