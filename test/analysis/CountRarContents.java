import java.io.*;
import java.net.HttpURLConnection;
import java.net.URL;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * Count files/folders directly under .rar files without downloading the entire archive.
 * Uses HTTP Range requests to download only the RAR header portion.
 */
public class CountRarContents {
	
	private static final int RAR_HEADER_SIZE = 1024 * 1024; // Download first 1MB (should contain file list)
	
	public static void main(String[] args) {
		if (args.length < 1) {
			System.out.println("Usage: java CountRarContents <rar_file_url>");
			System.out.println("Example: java CountRarContents https://ftp.pride.ebi.ac.uk/pride/data/archive/2016/08/PXD003791/file.rar");
			return;
		}
		
		String rarUrl = args[0];
		System.out.println("Analyzing RAR file: " + rarUrl);
		System.out.println("Downloading header only (first " + (RAR_HEADER_SIZE / 1024) + " KB)...");
		
		try {
			int count = countTopLevelEntriesInRar(rarUrl);
			System.out.println("\nTop-level entries in RAR file: " + count);
		} catch (Exception e) {
			System.err.println("Error: " + e.getMessage());
			e.printStackTrace();
		}
	}
	
	/**
	 * Count top-level files/folders in a RAR archive by reading only the header
	 */
	public static int countTopLevelEntriesInRar(String rarUrl) throws IOException {
		// Download RAR header using HTTP Range request
		byte[] headerData = downloadRarHeader(rarUrl);
		
		// Parse RAR header to extract file list
		Set<String> topLevelEntries = parseRarHeader(headerData);
		
		return topLevelEntries.size();
	}
	
	/**
	 * Download only the header portion of a RAR file using HTTP Range request
	 */
	private static byte[] downloadRarHeader(String urlString) throws IOException {
		URL url = new URL(urlString);
		HttpURLConnection connection = (HttpURLConnection) url.openConnection();
		
		// Use Range request to download only the first portion
		connection.setRequestProperty("Range", "bytes=0-" + (RAR_HEADER_SIZE - 1));
		connection.setRequestProperty("User-Agent", "Mozilla/5.0");
		connection.setConnectTimeout(15000);
		connection.setReadTimeout(15000);
		
		int responseCode = connection.getResponseCode();
		
		// 206 = Partial Content (expected for Range requests)
		// 200 = OK (server doesn't support Range, will get full file - but we'll only read first part)
		if (responseCode != HttpURLConnection.HTTP_PARTIAL && responseCode != HttpURLConnection.HTTP_OK) {
			throw new IOException("HTTP error: " + responseCode);
		}
		
		InputStream inputStream = connection.getInputStream();
		ByteArrayOutputStream buffer = new ByteArrayOutputStream();
		
		byte[] data = new byte[8192];
		int bytesRead;
		int totalRead = 0;
		
		while ((bytesRead = inputStream.read(data)) != -1 && totalRead < RAR_HEADER_SIZE) {
			buffer.write(data, 0, bytesRead);
			totalRead += bytesRead;
		}
		
		inputStream.close();
		connection.disconnect();
		
		System.out.println("Downloaded " + totalRead + " bytes of header");
		return buffer.toByteArray();
	}
	
	/**
	 * Parse RAR header to extract top-level file/folder names
	 * RAR format: Each file entry has a header block with the filename
	 */
	private static Set<String> parseRarHeader(byte[] headerData) {
		Set<String> topLevelEntries = new HashSet<>();
		
		try {
			// RAR file signature: Rar! (52 61 72 21 1A 07 00)
			// After signature, there are file header blocks
			
			// Simple approach: Look for file headers and extract filenames
			// RAR file header structure (simplified):
			// - Header CRC (2 bytes)
			// - Header type (1 byte)
			// - Header flags (2 bytes)
			// - Header size (2 bytes)
			// - Additional header data (variable)
			// - File name (null-terminated string)
			
			int pos = 0;
			
			// Find RAR signature
			boolean foundSignature = false;
			for (int i = 0; i < headerData.length - 7; i++) {
				if (headerData[i] == 0x52 && headerData[i+1] == 0x61 && headerData[i+2] == 0x72 && 
				    headerData[i+3] == 0x21 && headerData[i+4] == 0x1A && headerData[i+5] == 0x07) {
					pos = i + 7;
					foundSignature = true;
					break;
				}
			}
			
			if (!foundSignature) {
				System.out.println("Warning: RAR signature not found, trying alternative parsing...");
				pos = 0;
			}
			
			// Parse file headers
			while (pos < headerData.length - 20) {
				// Check if this looks like a file header
				// Header type: 0x74 = file header, 0x7A = archive header, 0x75 = subblock
				byte headerType = headerData[pos + 2];
				
				if (headerType == 0x74 || headerType == 0x7A || headerType == 0x75) {
					// Read header size (little-endian, 2 bytes at offset 5)
					if (pos + 7 >= headerData.length) break;
					
					int headerSize = (headerData[pos + 5] & 0xFF) | ((headerData[pos + 6] & 0xFF) << 8);
					
					if (headerSize < 7 || headerSize > 65535 || pos + headerSize > headerData.length) {
						pos++;
						continue;
					}
					
					// Try to extract filename (usually starts around offset 7-20 in header)
					// Filename is null-terminated string
					for (int nameStart = pos + 7; nameStart < pos + headerSize && nameStart < headerData.length - 1; nameStart++) {
						// Look for printable ASCII characters (likely filename start)
						if (headerData[nameStart] >= 32 && headerData[nameStart] < 127) {
							StringBuilder filename = new StringBuilder();
							int namePos = nameStart;
							
							// Read until null terminator or non-printable
							while (namePos < pos + headerSize && namePos < headerData.length) {
								byte b = headerData[namePos];
								if (b == 0) break; // Null terminator
								if (b < 32 || b >= 127) break; // Non-printable
								filename.append((char) b);
								namePos++;
							}
							
							if (filename.length() > 0 && filename.length() < 256) {
								String name = filename.toString();
								// Extract top-level name (before first / or \)
								if (name.contains("/")) {
									name = name.substring(0, name.indexOf('/'));
								} else if (name.contains("\\")) {
									name = name.substring(0, name.indexOf('\\'));
								}
								
								// Only count if it looks like a valid filename
								if (name.length() > 0 && !name.equals(".") && !name.equals("..")) {
									topLevelEntries.add(name);
								}
							}
							
							// Move past this filename
							pos = namePos;
							break;
						}
					}
					
					pos += headerSize;
				} else {
					pos++;
				}
			}
			
		} catch (Exception e) {
			System.err.println("Error parsing RAR header: " + e.getMessage());
			// Fallback: try to find filenames by looking for common patterns
			topLevelEntries.addAll(findFilenamesByPattern(headerData));
		}
		
		return topLevelEntries;
	}
	
	/**
	 * Fallback method: Try to find filenames by looking for patterns
	 */
	private static Set<String> findFilenamesByPattern(byte[] data) {
		Set<String> entries = new HashSet<>();
		
		// Look for sequences of printable ASCII that might be filenames
		StringBuilder current = new StringBuilder();
		
		for (int i = 0; i < data.length; i++) {
			byte b = data[i];
			if (b >= 32 && b < 127 && b != 0) {
				current.append((char) b);
			} else {
				if (current.length() > 3 && current.length() < 256) {
					String str = current.toString();
					// Check if it looks like a filename (contains extension or common patterns)
					if (str.contains(".") || str.matches("^[A-Za-z0-9_\\-]+$")) {
						// Extract top-level name
						if (str.contains("/")) {
							str = str.substring(0, str.indexOf('/'));
						} else if (str.contains("\\")) {
							str = str.substring(0, str.indexOf('\\'));
						}
						if (str.length() > 0) {
							entries.add(str);
						}
					}
				}
				current.setLength(0);
			}
		}
		
		return entries;
	}
}

