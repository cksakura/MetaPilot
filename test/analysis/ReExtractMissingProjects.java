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
 * Re-extract .raw files for missing projects and append to existing raw_files_list.txt
 */
public class ReExtractMissingProjects {
	
	public static void main(String[] args) {
		// Try multiple possible paths - use absolute paths to avoid duplication
		File missingProjectsFile = null;
		File existingRawFilesFile = null;
		File outputFile = null;
		
		String[] possiblePaths = {
			"missing_projects.txt",
			"./missing_projects.txt",
			"test/analysis/missing_projects.txt",
			"analysis/missing_projects.txt"
		};
		
		for (String path : possiblePaths) {
			File f = new File(path);
			if (f.exists()) {
				missingProjectsFile = f.getAbsoluteFile();
				// Use the same directory as missing_projects.txt for other files
				File baseDir = missingProjectsFile.getParentFile();
				if (baseDir != null) {
					existingRawFilesFile = new File(baseDir, "raw_files_list.txt").getAbsoluteFile();
					outputFile = new File(baseDir, "raw_files_list_complete.txt").getAbsoluteFile();
				} else {
					// If no parent, use current directory
					existingRawFilesFile = new File("raw_files_list.txt").getAbsoluteFile();
					outputFile = new File("raw_files_list_complete.txt").getAbsoluteFile();
				}
				break;
			}
		}
		
		if (missingProjectsFile == null || !missingProjectsFile.exists()) {
			System.err.println("Error: missing_projects.txt not found");
			System.err.println("Searched in current directory and common paths");
			return;
		}
		
		System.out.println("Reading missing projects from: " + missingProjectsFile.getAbsolutePath());
		System.out.println("Reading existing raw files from: " + existingRawFilesFile.getAbsolutePath());
		System.out.println("Writing complete list to: " + outputFile.getAbsolutePath());
		
		try {
			// Read existing raw files
			Set<String> existingEntries = new HashSet<>();
			if (existingRawFilesFile.exists()) {
				BufferedReader reader = new BufferedReader(new FileReader(existingRawFilesFile));
				String line = reader.readLine(); // Skip header
				while ((line = reader.readLine()) != null) {
					line = line.trim();
					if (!line.isEmpty()) {
						existingEntries.add(line);
					}
				}
				reader.close();
				System.out.println("Found " + existingEntries.size() + " existing entries");
			}
			
			// Read missing projects
			List<String> missingProjectLinks = new ArrayList<>();
			BufferedReader reader = new BufferedReader(new FileReader(missingProjectsFile));
			String line = reader.readLine(); // Skip header
			line = reader.readLine(); // Skip separator
			line = reader.readLine(); // Skip column header
			line = reader.readLine(); // Skip separator
			
			while ((line = reader.readLine()) != null) {
				line = line.trim();
				if (line.isEmpty() || line.startsWith("=")) {
					continue;
				}
				
				String[] parts = line.split("\t");
				if (parts.length >= 2) {
					String ftpLink = parts[1].trim();
					if (!ftpLink.isEmpty()) {
						missingProjectLinks.add(ftpLink);
					}
				}
			}
			reader.close();
			
			System.out.println("Found " + missingProjectLinks.size() + " missing projects to re-extract");
			
			// Re-extract .raw files for missing projects
			List<String> newRawFiles = new ArrayList<>();
			
			for (int i = 0; i < missingProjectLinks.size(); i++) {
				String ftpLink = missingProjectLinks.get(i);
				System.out.println("\nRe-extracting " + (i + 1) + "/" + missingProjectLinks.size() + ": " + ftpLink);
				
				List<String> rawFiles = new ArrayList<>();
				
				// Try HTTPS first (original URL)
				rawFiles = extractRawFilesFromHttp(ftpLink, true);
				
				// If HTTPS fails, try HTTP
				if (rawFiles.isEmpty()) {
					System.out.println("  HTTPS failed, trying HTTP...");
					rawFiles = extractRawFilesFromHttp(ftpLink, false);
				}
				
				// If HTTP/HTTPS fails, try FTP
				if (rawFiles.isEmpty()) {
					System.out.println("  HTTP/HTTPS failed, trying FTP...");
					rawFiles = extractRawFilesFromFtp(ftpLink);
				}
				
				if (rawFiles.isEmpty()) {
					System.out.println("  WARNING: Still no .raw files found!");
				} else {
					System.out.println("  Found " + rawFiles.size() + " .raw file(s)");
					for (String rawFile : rawFiles) {
						String entry = ftpLink + "\t" + rawFile;
						if (!existingEntries.contains(entry)) {
							newRawFiles.add(entry);
						}
					}
				}
				
				// Small delay to avoid overwhelming the server
				Thread.sleep(500);
			}
			
			// Ensure output directory exists
			File outputDir = outputFile.getParentFile();
			if (outputDir != null && !outputDir.exists()) {
				outputDir.mkdirs();
			}
			
			// Write complete list
			PrintWriter writer = new PrintWriter(outputFile);
			writer.println("FTP_Link\tRaw_File_Name");
			
			// Write existing entries
			for (String entry : existingEntries) {
				writer.println(entry);
			}
			
			// Write new entries
			for (String entry : newRawFiles) {
				writer.println(entry);
			}
			
			writer.close();
			
			System.out.println("\n" + "=".repeat(60));
			System.out.println("Re-extraction complete!");
			System.out.println("New .raw files found: " + newRawFiles.size());
			System.out.println("Total entries in complete list: " + (existingEntries.size() + newRawFiles.size()));
			System.out.println("Complete list written to: " + outputFile.getAbsolutePath());
			
		} catch (Exception e) {
			System.err.println("Error: " + e.getMessage());
			e.printStackTrace();
		}
	}
	
	// Copy the HTTP extraction method from ExtractRawFiles
	private static List<String> extractRawFilesFromHttp(String urlString, boolean useHttps) {
		List<String> rawFiles = new ArrayList<>();
		
		try {
			String urlToUse = urlString;
			if (!useHttps) {
				urlToUse = urlString.replace("https://", "http://");
			}
			
			java.net.URL url = new java.net.URL(urlToUse);
			java.net.URLConnection connection = url.openConnection();
			connection.setConnectTimeout(15000);
			connection.setReadTimeout(15000);
			connection.setRequestProperty("User-Agent", "Mozilla/5.0");
			
			BufferedReader reader = new BufferedReader(
				new java.io.InputStreamReader(connection.getInputStream()));
			
			String line;
			Pattern rawPattern1 = Pattern.compile("href\\s*=\\s*[\"']([^\"']+\\.raw)[\"']", Pattern.CASE_INSENSITIVE);
			Pattern rawPattern2 = Pattern.compile("href\\s*=\\s*[\"']([^\"']+\\.raw)", Pattern.CASE_INSENSITIVE);
			Pattern rawPattern3 = Pattern.compile(">\\s*([^<>\\s]+\\.raw)\\s*<", Pattern.CASE_INSENSITIVE);
			Pattern rawPattern4 = Pattern.compile("\\b([^\\s<>\"']+\\.raw)\\b", Pattern.CASE_INSENSITIVE);
			
			StringBuilder content = new StringBuilder();
			while ((line = reader.readLine()) != null) {
				content.append(line).append("\n");
				
				Matcher matcher = rawPattern1.matcher(line);
				while (matcher.find()) {
					String rawFile = matcher.group(1);
					if (!rawFiles.contains(rawFile)) {
						rawFiles.add(rawFile);
					}
				}
				
				matcher = rawPattern2.matcher(line);
				while (matcher.find()) {
					String rawFile = matcher.group(1);
					if (!rawFiles.contains(rawFile)) {
						rawFiles.add(rawFile);
					}
				}
				
				matcher = rawPattern3.matcher(line);
				while (matcher.find()) {
					String rawFile = matcher.group(1);
					if (!rawFiles.contains(rawFile)) {
						rawFiles.add(rawFile);
					}
				}
			}
			reader.close();
			
			if (rawFiles.isEmpty()) {
				Pattern ftpPattern = Pattern.compile("\\S+\\s+\\d+\\s+\\S+\\s+\\S+\\s+\\d+\\s+\\S+\\s+\\d+\\s+\\S+\\s+(.+\\.raw)", Pattern.CASE_INSENSITIVE);
				Matcher matcher = ftpPattern.matcher(content.toString());
				while (matcher.find()) {
					String rawFile = matcher.group(1).trim();
					if (!rawFiles.contains(rawFile)) {
						rawFiles.add(rawFile);
					}
				}
				
				if (rawFiles.isEmpty()) {
					matcher = rawPattern4.matcher(content.toString());
					while (matcher.find()) {
						String rawFile = matcher.group(1).trim();
						if (rawFile.contains("/") || rawFile.contains("\\")) {
							rawFile = rawFile.substring(rawFile.lastIndexOf('/') + 1);
							rawFile = rawFile.substring(rawFile.lastIndexOf('\\') + 1);
						}
						if (rawFile.endsWith(".raw") && !rawFiles.contains(rawFile)) {
							rawFiles.add(rawFile);
						}
					}
				}
			}
			
		} catch (Exception e) {
			// Silent failure, will try other methods
		}
		
		return rawFiles;
	}
	
	// Copy the FTP extraction method from ExtractRawFiles
	private static List<String> extractRawFilesFromFtp(String urlString) {
		List<String> rawFiles = new ArrayList<>();
		
		try {
			Class<?> ftpClientClass = Class.forName("org.apache.commons.net.ftp.FTPClient");
			Class<?> ftpFileClass = Class.forName("org.apache.commons.net.ftp.FTPFile");
			
			Object ftpClient = ftpClientClass.getDeclaredConstructor().newInstance();
			
			java.net.URL url = new java.net.URL(urlString);
			String host = url.getHost();
			String path = url.getPath();
			
			java.lang.reflect.Method connectMethod = ftpClientClass.getMethod("connect", String.class);
			connectMethod.invoke(ftpClient, host);
			
			java.lang.reflect.Method enterPassiveMethod = ftpClientClass.getMethod("enterLocalPassiveMode");
			enterPassiveMethod.invoke(ftpClient);
			
			java.lang.reflect.Method loginMethod = ftpClientClass.getMethod("login", String.class, String.class);
			loginMethod.invoke(ftpClient, "anonymous", "");
			
			if (path != null && !path.isEmpty()) {
				java.lang.reflect.Method changeDirMethod = ftpClientClass.getMethod("changeWorkingDirectory", String.class);
				changeDirMethod.invoke(ftpClient, path);
			}
			
			java.lang.reflect.Method listFilesMethod = ftpClientClass.getMethod("listFiles");
			Object[] files = (Object[]) listFilesMethod.invoke(ftpClient);
			
			java.lang.reflect.Method isFileMethod = ftpFileClass.getMethod("isFile");
			java.lang.reflect.Method getNameMethod = ftpFileClass.getMethod("getName");
			
			for (Object file : files) {
				Boolean isFile = (Boolean) isFileMethod.invoke(file);
				if (isFile) {
					String fileName = (String) getNameMethod.invoke(file);
					if (fileName.toLowerCase().endsWith(".raw")) {
						rawFiles.add(fileName);
					}
				}
			}
			
			java.lang.reflect.Method isConnectedMethod = ftpClientClass.getMethod("isConnected");
			Boolean isConnected = (Boolean) isConnectedMethod.invoke(ftpClient);
			if (isConnected) {
				java.lang.reflect.Method disconnectMethod = ftpClientClass.getMethod("disconnect");
				disconnectMethod.invoke(ftpClient);
			}
			
		} catch (ClassNotFoundException e) {
			System.out.println("  FTP not available (Apache Commons Net not found)");
		} catch (Exception e) {
			// FTP failed
		}
		
		return rawFiles;
	}
	
	private static String repeatString(String str, int count) {
		StringBuilder sb = new StringBuilder();
		for (int i = 0; i < count; i++) {
			sb.append(str);
		}
		return sb.toString();
	}
}

