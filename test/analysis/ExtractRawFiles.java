import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.net.URL;
import java.net.URLConnection;
import java.util.ArrayList;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
// FTP support (optional - will gracefully degrade if not available)
// import org.apache.commons.net.ftp.FTPClient;
// import org.apache.commons.net.ftp.FTPFile;

/**
 * Extract .raw file names from PRIDE FTP/HTTP directory links
 */
public class ExtractRawFiles {
	
	public static void main(String[] args) {
		// Try multiple possible paths for projects.txt
		File projectsFile = null;
		String[] possiblePaths = {
			"test/analysis/projects.txt",
			"projects.txt",
			"./projects.txt",
			"analysis/projects.txt"
		};
		
		for (String path : possiblePaths) {
			File f = new File(path);
			if (f.exists()) {
				projectsFile = f;
				break;
			}
		}
		
		if (projectsFile == null) {
			System.err.println("Error: projects.txt not found. Tried:");
			for (String path : possiblePaths) {
				System.err.println("  " + new File(path).getAbsolutePath());
			}
			return;
		}
		
		File outputFile = new File(projectsFile.getParent(), "raw_files_list.txt");
		System.out.println("Reading from: " + projectsFile.getAbsolutePath());
		System.out.println("Writing to: " + outputFile.getAbsolutePath());
		
		try {
			List<String> ftpLinks = readFtpLinks(projectsFile);
			System.out.println("Found " + ftpLinks.size() + " FTP links to process");
			
			List<String> allRawFiles = new ArrayList<>();
			
			for (int i = 0; i < ftpLinks.size(); i++) {
				String ftpLink = ftpLinks.get(i);
				System.out.println("Processing " + (i + 1) + "/" + ftpLinks.size() + ": " + ftpLink);
				
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
				
				// Check for .rar files that might contain .raw files
				List<String> rarFiles = findRarFiles(ftpLink);
				if (!rarFiles.isEmpty()) {
					System.out.println("  Found " + rarFiles.size() + " .rar file(s)");
					for (String rarFile : rarFiles) {
						String rarUrl = ftpLink + (ftpLink.endsWith("/") ? "" : "/") + rarFile;
						try {
							int entryCount = CountRarContents.countTopLevelEntriesInRar(rarUrl);
							System.out.println("    " + rarFile + " contains " + entryCount + " top-level entries");
							// Note: We don't extract individual .raw files from RAR, just report the count
							allRawFiles.add(ftpLink + "\t" + rarFile + " (RAR with " + entryCount + " entries)");
						} catch (Exception e) {
							System.out.println("    " + rarFile + " - Could not read RAR contents: " + e.getMessage());
							allRawFiles.add(ftpLink + "\t" + rarFile + " (RAR - unable to read)");
						}
					}
				}
				
				if (rawFiles.isEmpty() && rarFiles.isEmpty()) {
					System.out.println("  WARNING: No .raw files found (may need manual check)");
				} else if (!rawFiles.isEmpty()) {
					System.out.println("  Found " + rawFiles.size() + " .raw file(s)");
					for (String rawFile : rawFiles) {
						allRawFiles.add(ftpLink + "\t" + rawFile);
					}
				}
				
				// Small delay to avoid overwhelming the server
				Thread.sleep(500);
			}
			
			// Write results
			java.io.PrintWriter writer = new java.io.PrintWriter(outputFile);
			writer.println("FTP_Link\tRaw_File_Name");
			for (String line : allRawFiles) {
				writer.println(line);
			}
			writer.close();
			
			System.out.println("\nResults written to: " + outputFile.getAbsolutePath());
			System.out.println("Total .raw files found: " + allRawFiles.size());
			
		} catch (Exception e) {
			System.err.println("Fatal error: " + e.getMessage());
			e.printStackTrace();
		}
	}
	
	/**
	 * Find .rar files in the directory listing
	 */
	private static List<String> findRarFiles(String urlString) {
		List<String> rarFiles = new ArrayList<>();
		
		try {
			// Try HTTPS first
			List<String> files = listFilesInDirectory(urlString, true);
			for (String file : files) {
				if (file.toLowerCase().endsWith(".rar")) {
					rarFiles.add(file);
				}
			}
			
			// If HTTPS fails, try HTTP
			if (rarFiles.isEmpty()) {
				files = listFilesInDirectory(urlString, false);
				for (String file : files) {
					if (file.toLowerCase().endsWith(".rar")) {
						rarFiles.add(file);
					}
				}
			}
		} catch (Exception e) {
			// Silent failure
		}
		
		return rarFiles;
	}
	
	/**
	 * List all files in a directory (helper method)
	 */
	private static List<String> listFilesInDirectory(String urlString, boolean useHttps) {
		List<String> files = new ArrayList<>();
		
		try {
			String urlToUse = urlString;
			if (!useHttps) {
				urlToUse = urlString.replace("https://", "http://");
			}
			
			URL url = new URL(urlToUse);
			URLConnection connection = url.openConnection();
			connection.setConnectTimeout(15000);
			connection.setReadTimeout(15000);
			connection.setRequestProperty("User-Agent", "Mozilla/5.0");
			
			BufferedReader reader = new BufferedReader(
				new InputStreamReader(connection.getInputStream()));
			
			String line;
			Pattern filePattern = Pattern.compile("href\\s*=\\s*[\"']([^\"']+)[\"']", Pattern.CASE_INSENSITIVE);
			
			while ((line = reader.readLine()) != null) {
				Matcher matcher = filePattern.matcher(line);
				while (matcher.find()) {
					String fileName = matcher.group(1);
					// Skip parent directory links
					if (!fileName.equals("../") && !fileName.equals("./") && !fileName.endsWith("/")) {
						files.add(fileName);
					}
				}
			}
			reader.close();
		} catch (Exception e) {
			// Silent failure
		}
		
		return files;
	}
	
	private static List<String> readFtpLinks(File file) throws IOException {
		List<String> links = new ArrayList<>();
		BufferedReader reader = new BufferedReader(new FileReader(file));
		String line;
		while ((line = reader.readLine()) != null) {
			line = line.trim();
			if (!line.isEmpty() && !line.startsWith("#")) {
				links.add(line);
			}
		}
		reader.close();
		return links;
	}
	
	/**
	 * Try to extract .raw files from HTTP/HTTPS URL
	 * @param useHttps If true, use HTTPS; if false, convert to HTTP
	 */
	private static List<String> extractRawFilesFromHttp(String urlString, boolean useHttps) {
		List<String> rawFiles = new ArrayList<>();
		
		try {
			// Use HTTPS or HTTP based on parameter
			String urlToUse = urlString;
			if (!useHttps) {
				urlToUse = urlString.replace("https://", "http://");
			}
			
			URL url = new URL(urlToUse);
			URLConnection connection = url.openConnection();
			connection.setConnectTimeout(15000); // 15 seconds
			connection.setReadTimeout(15000);
			connection.setRequestProperty("User-Agent", "Mozilla/5.0");
			
			// Try to read as HTML directory listing
			BufferedReader reader = new BufferedReader(
				new InputStreamReader(connection.getInputStream()));
			
			String line;
			// More comprehensive patterns to match .raw file links in HTML
			Pattern rawPattern1 = Pattern.compile("href\\s*=\\s*[\"']([^\"']+\\.raw)[\"']", Pattern.CASE_INSENSITIVE);
			Pattern rawPattern2 = Pattern.compile("href\\s*=\\s*[\"']([^\"']+\\.raw)", Pattern.CASE_INSENSITIVE);
			Pattern rawPattern3 = Pattern.compile(">\\s*([^<>\\s]+\\.raw)\\s*<", Pattern.CASE_INSENSITIVE);
			Pattern rawPattern4 = Pattern.compile("\\b([^\\s<>\"']+\\.raw)\\b", Pattern.CASE_INSENSITIVE);
			
			StringBuilder content = new StringBuilder();
			while ((line = reader.readLine()) != null) {
				content.append(line).append("\n");
				
				// Try pattern 1: href="filename.raw"
				Matcher matcher = rawPattern1.matcher(line);
				while (matcher.find()) {
					String rawFile = matcher.group(1);
					if (!rawFiles.contains(rawFile)) {
						rawFiles.add(rawFile);
					}
				}
				
				// Try pattern 2: href="filename.raw (without closing quote)
				matcher = rawPattern2.matcher(line);
				while (matcher.find()) {
					String rawFile = matcher.group(1);
					if (!rawFiles.contains(rawFile)) {
						rawFiles.add(rawFile);
					}
				}
				
				// Try pattern 3: >filename.raw<
				matcher = rawPattern3.matcher(line);
				while (matcher.find()) {
					String rawFile = matcher.group(1);
					if (!rawFiles.contains(rawFile)) {
						rawFiles.add(rawFile);
					}
				}
			}
			reader.close();
			
			// If no matches found with HTML patterns, try parsing FTP-style listing or plain text
			if (rawFiles.isEmpty()) {
				// FTP LIST format: permissions size date filename
				Pattern ftpPattern = Pattern.compile("\\S+\\s+\\d+\\s+\\S+\\s+\\S+\\s+\\d+\\s+\\S+\\s+\\d+\\s+\\S+\\s+(.+\\.raw)", Pattern.CASE_INSENSITIVE);
				Matcher matcher = ftpPattern.matcher(content.toString());
				while (matcher.find()) {
					String rawFile = matcher.group(1).trim();
					if (!rawFiles.contains(rawFile)) {
						rawFiles.add(rawFile);
					}
				}
				
				// If still empty, try a more general pattern for any .raw file mention
				if (rawFiles.isEmpty()) {
					matcher = rawPattern4.matcher(content.toString());
					while (matcher.find()) {
						String rawFile = matcher.group(1).trim();
						// Filter out URLs and paths, just get filenames
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
			// HTTP failed, will try other methods
			// Uncomment for debugging: System.err.println("  HTTP failed: " + e.getMessage());
		}
		
		return rawFiles;
	}
	
	/**
	 * Try to extract .raw files from FTP URL
	 * Note: Requires Apache Commons Net library
	 */
	private static List<String> extractRawFilesFromFtp(String urlString) {
		List<String> rawFiles = new ArrayList<>();
		
		try {
			// Use reflection to check if FTP classes are available
			Class<?> ftpClientClass = Class.forName("org.apache.commons.net.ftp.FTPClient");
			Class<?> ftpFileClass = Class.forName("org.apache.commons.net.ftp.FTPFile");
			
			Object ftpClient = ftpClientClass.getDeclaredConstructor().newInstance();
			
			// Parse URL to get host and path
			URL url = new URL(urlString);
			String host = url.getHost();
			String path = url.getPath();
			
			// Connect to FTP server using reflection
			java.lang.reflect.Method connectMethod = ftpClientClass.getMethod("connect", String.class);
			connectMethod.invoke(ftpClient, host);
			
			java.lang.reflect.Method enterPassiveMethod = ftpClientClass.getMethod("enterLocalPassiveMode");
			enterPassiveMethod.invoke(ftpClient);
			
			java.lang.reflect.Method loginMethod = ftpClientClass.getMethod("login", String.class, String.class);
			loginMethod.invoke(ftpClient, "anonymous", "");
			
			// Change to the directory
			if (path != null && !path.isEmpty()) {
				java.lang.reflect.Method changeDirMethod = ftpClientClass.getMethod("changeWorkingDirectory", String.class);
				changeDirMethod.invoke(ftpClient, path);
			}
			
			// List files
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
			
			// Disconnect
			java.lang.reflect.Method isConnectedMethod = ftpClientClass.getMethod("isConnected");
			Boolean isConnected = (Boolean) isConnectedMethod.invoke(ftpClient);
			if (isConnected) {
				java.lang.reflect.Method disconnectMethod = ftpClientClass.getMethod("disconnect");
				disconnectMethod.invoke(ftpClient);
			}
			
		} catch (ClassNotFoundException e) {
			// Apache Commons Net not available - FTP disabled
			System.out.println("  FTP not available (Apache Commons Net not found), skipping FTP method");
		} catch (Exception e) {
			// FTP failed, return empty list
			// Uncomment for debugging: System.err.println("  FTP failed: " + e.getMessage());
		}
		
		return rawFiles;
	}
}

