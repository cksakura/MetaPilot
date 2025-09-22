/**
 * 
 */
package uk.ac.quadram.metapilot.mdb.unipept;

import java.io.BufferedReader;
import java.io.DataOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.net.HttpURLConnection;
import java.net.MalformedURLException;
import java.net.ProtocolException;
import java.net.URL;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

/**
 * @author Kai Cheng
 *
 */
public class UnipeptRequester {

	private static final String USER_AGENT = "Mozilla/5.0 (Windows NT 6.1; WOW64) ...";
	private static final String Language = "en-US,en;q=0.5";
	private static final String pept2prot = "http://api.unipept.ugent.be/api/v2/pept2prot";
	private static final String pept2lca = "http://api.unipept.ugent.be/api/v2/pept2lca";
	private static final String protinfo = "http://api.unipept.ugent.be/api/v2/protinfo";
	// private static final String pept2taxa =
	// "http://api.unipept.ugent.be/api/v1/pept2taxa";

	private static Logger LOGGER = LogManager.getLogger(UnipeptRequester.class);

	public static UnipProResult pept2prot(String sequence) {

		URL obj = null;
		try {
			obj = new URL(pept2prot);
		} catch (MalformedURLException e) {
			// TODO Auto-generated catch block
			LOGGER.error("Error in matching peptide " + sequence + " to proteins by " + pept2prot, e);
			return null;
		}
		HttpURLConnection con = null;
		try {
			con = (HttpURLConnection) obj.openConnection();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			LOGGER.error("Error in matching peptide " + sequence + " to proteins by " + pept2prot, e);
			return null;
		}

		// add reuqest header
		try {
			con.setRequestMethod("POST");
		} catch (ProtocolException e) {
			// TODO Auto-generated catch block
			LOGGER.error("Error in matching peptide " + sequence + " to proteins by " + pept2prot, e);
		}

		con.setRequestProperty("accept", "*/*");
		con.setRequestProperty("connection", "Keep-Alive");
		con.setRequestProperty("User-Agent", USER_AGENT);
		con.setRequestProperty("Accept-Language", Language);

		String urlParameters = "input[]=" + sequence;

		// Send post request
		con.setDoOutput(true);
		con.setDoInput(true);

		DataOutputStream wr;
		try {
			wr = new DataOutputStream(con.getOutputStream());
			wr.writeBytes(urlParameters);
			wr.flush();
			wr.close();

		} catch (IOException e) {
			// TODO Auto-generated catch block
			LOGGER.error("Error in matching peptide " + sequence + " to proteins by " + pept2prot, e);
		}

		// int responseCode = con.getResponseCode();
		// System.out.println("\nSending 'POST' request to URL : " + pept2prot);
		// System.out.println("Post parameters : " + urlParameters);
		// System.out.println("Response Code : " + responseCode);

		String line = null;
		try {
			BufferedReader in = new BufferedReader(new InputStreamReader(con.getInputStream()));
			line = in.readLine();
			in.close();

		} catch (IOException e) {
			// TODO Auto-generated catch block
			LOGGER.error("Error in matching peptide " + sequence + " to proteins by " + pept2prot, e);
		}

		if (line == null || line.length() <= 2) {
			return null;
		}

		try {
			// Remove the outer brackets and split by },{
			String jsonContent = line.substring(1, line.length() - 1);
			String[] entries = jsonContent.split("\\},\\{");
			
			String[] uniprotIds = new String[entries.length];
			String[] pronames = new String[entries.length];
			int[] taxonIds = new int[entries.length];

			for (int i = 0; i < entries.length; i++) {
				String entry = entries[i];
				
				// Parse JSON-like structure more carefully
				// Expected format: {"peptide":"...","uniprot_id":"...","protein_name":"...","taxon_id":...}
				uniprotIds[i] = extractJsonValue(entry, "uniprot_id");
				pronames[i] = extractJsonValue(entry, "protein_name");
				
				String taxonIdStr = extractJsonValue(entry, "taxon_id");
				try {
					taxonIds[i] = Integer.parseInt(taxonIdStr);
				} catch (NumberFormatException e) {
					LOGGER.warn("Error parsing taxon ID '" + taxonIdStr + "': " + e.getMessage());
					taxonIds[i] = 0; // Default value
				}
			}

			UnipProResult result = new UnipProResult(sequence, uniprotIds, pronames, taxonIds);
			return result;
			
		} catch (Exception e) {
			LOGGER.error("Error in parsing response: " + e.getMessage(), e);
			return null;
		}
	}

	public static HashMap<String, UnipProResult> pept2prot(String[] sequences) {

		URL obj = null;
		try {
			obj = new URL(pept2prot);
		} catch (MalformedURLException e) {
			// TODO Auto-generated catch block
			LOGGER.error("Error in matching peptides to proteins by " + pept2prot, e);
			return null;
		}
		HttpURLConnection con = null;
		try {
			con = (HttpURLConnection) obj.openConnection();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			LOGGER.error("Error in matching peptides to proteins by " + pept2prot, e);
			return null;
		}

		// add reuqest header
		try {
			con.setRequestMethod("POST");
		} catch (ProtocolException e) {
			// TODO Auto-generated catch block
			LOGGER.error("Error in matching peptides to proteins by " + pept2prot, e);
		}
		con.setRequestProperty("accept", "*/*");
		con.setRequestProperty("connection", "Keep-Alive");
		con.setRequestProperty("User-Agent", USER_AGENT);
		con.setRequestProperty("Accept-Language", Language);

		StringBuilder sb = new StringBuilder();
		for (String sequence : sequences) {
			sb.append("input[]=").append(sequence).append("&");
		}
		sb.deleteCharAt(sb.length() - 1);

		String urlParameters = sb.toString();

		// Send post request
		con.setDoOutput(true);
		con.setDoInput(true);

		try {
			DataOutputStream wr = new DataOutputStream(con.getOutputStream());
			wr.writeBytes(urlParameters);
			wr.flush();
			wr.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			LOGGER.error("Error in matching peptides to proteins by " + pept2prot, e);
		}

		HashMap<String, UnipProResult> map = new HashMap<String, UnipProResult>();

		String line = null;
		try {
			BufferedReader in = new BufferedReader(new InputStreamReader(con.getInputStream()));
			line = in.readLine();
			in.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			LOGGER.error("Error in matching peptides to proteins by " + pept2prot, e);
		}

		if (line == null || line.length() < 4) {
			return map;
		}

		HashMap<String, HashSet<String>> pepProMap = new HashMap<String, HashSet<String>>();
		HashMap<String, String[]> proTaxonMap = new HashMap<String, String[]>();

		try {
			// Remove the outer brackets and split by },{
			String jsonContent = line.substring(1, line.length() - 1);
			String[] entries = jsonContent.split("\\},\\{");

			for (int i = 0; i < entries.length; i++) {
				String entry = entries[i];
				
				// Parse JSON-like structure more carefully
				// Expected format: {"peptide":"...","uniprot_id":"...","protein_name":"...","taxon_id":...}
				String sequence = extractJsonValue(entry, "peptide");
				String proid = extractJsonValue(entry, "uniprot_id");
				String proname = extractJsonValue(entry, "protein_name");
				String taxid = extractJsonValue(entry, "taxon_id");

				if (!sequence.isEmpty() && !proid.isEmpty()) {
					if (pepProMap.containsKey(sequence)) {
						pepProMap.get(sequence).add(proid);
					} else {
						HashSet<String> set = new HashSet<String>();
						set.add(proid);
						pepProMap.put(sequence, set);
					}

					proTaxonMap.put(proid, new String[] { proname, taxid });
				}
			}

			for (String sequence : pepProMap.keySet()) {
				HashSet<String> set = pepProMap.get(sequence);
				String[] proids = set.toArray(new String[set.size()]);
				String[] pronames = new String[proids.length];
				int[] taxids = new int[proids.length];

				for (int i = 0; i < proids.length; i++) {
					String[] proTaxonInfo = proTaxonMap.get(proids[i]);
					if (proTaxonInfo != null && proTaxonInfo.length >= 2) {
						pronames[i] = proTaxonInfo[0];
						try {
							taxids[i] = Integer.parseInt(proTaxonInfo[1]);
						} catch (NumberFormatException e) {
							LOGGER.warn("Error parsing taxon ID '" + proTaxonInfo[1] + "': " + e.getMessage());
							taxids[i] = 0; // Default value
						}
					} else {
						pronames[i] = "";
						taxids[i] = 0;
					}
				}

				UnipProResult result = new UnipProResult(sequence, proids, pronames, taxids);
				map.put(sequence, result);
			}
			
		} catch (Exception e) {
			LOGGER.error("Error in parsing response: " + e.getMessage(), e);
		}

		return map;
	}

	public static UnipLCAResult pept2lca(String sequence) {
		LOGGER.debug("Starting pept2lca request for single sequence: " + sequence);

		URL obj = null;
		try {
			obj = new URL(pept2lca);
		} catch (MalformedURLException e) {
			// TODO Auto-generated catch block
			LOGGER.error("Error in matching peptide " + sequence + " to LCA by " + pept2lca, e);
			return null;
		}
		HttpURLConnection con = null;
		try {
			con = (HttpURLConnection) obj.openConnection();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			LOGGER.error("Error in matching peptide " + sequence + " to LCA by " + pept2lca, e);
			return null;
		}

		// add reuqest header
		try {
			con.setRequestMethod("POST");
		} catch (ProtocolException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		con.setRequestProperty("accept", "*/*");
		con.setRequestProperty("connection", "Keep-Alive");
		con.setRequestProperty("User-Agent", USER_AGENT);
		con.setRequestProperty("Accept-Language", Language);

		String urlParameters = "input[]=" + sequence;

		// Send post request
		con.setDoOutput(true);
		con.setDoInput(true);

		try {
			DataOutputStream wr = new DataOutputStream(con.getOutputStream());
			wr.writeBytes(urlParameters);
			wr.flush();
			wr.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			LOGGER.error("Error in matching peptide " + sequence + " to LCA by " + pept2lca, e);
		}

		// int responseCode = con.getResponseCode();
		// System.out.println("\nSending 'POST' request to URL : " + pept2prot);
		// System.out.println("Post parameters : " + urlParameters);
		// System.out.println("Response Code : " + responseCode);

		String line = null;
		try {
			BufferedReader in = new BufferedReader(new InputStreamReader(con.getInputStream()));
			line = in.readLine();
			in.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			LOGGER.error("Error in matching peptide " + sequence + " to LCA by " + pept2lca, e);
		}

		if (line == null || line.length() == 0) {
			LOGGER.warn("Empty response from Unipept API for sequence: " + sequence);
			return null;
		}

		String[] cs = line.split("[\\W]+");
		if (cs.length != 9) {
			LOGGER.warn("Unexpected response format for sequence " + sequence + ": expected 9 fields, got " + cs.length + ". Response: " + line);
			return null;
		}

		UnipLCAResult result = new UnipLCAResult(cs[2], Integer.parseInt(cs[4]), cs[6], cs[8]);
		LOGGER.debug("Successfully parsed LCA result for sequence " + sequence + ": taxon ID " + result.getTaxon().getId());
		return result;
	}

	public static HashMap<String, UnipLCAResult> pept2lca(String[] sequences) {
		LOGGER.info("Starting pept2lca request for " + sequences.length + " sequences");

		URL obj = null;
		try {
			obj = new URL(pept2lca);
		} catch (MalformedURLException e) {
			// TODO Auto-generated catch block
			LOGGER.error("Error in matching peptides to LCA by " + pept2lca, e);
			return null;
		}
		HttpURLConnection con = null;
		try {
			con = (HttpURLConnection) obj.openConnection();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			LOGGER.error("Error in matching peptides to LCA by " + pept2lca, e);
			return null;
		}

		// add reuqest header
		try {
			con.setRequestMethod("POST");
		} catch (ProtocolException e) {
			// TODO Auto-generated catch block
			LOGGER.error("Error in matching peptides to LCA by " + pept2lca, e);
		}
		con.setRequestProperty("accept", "*/*");
		con.setRequestProperty("connection", "Keep-Alive");
		con.setRequestProperty("User-Agent", USER_AGENT);
		con.setRequestProperty("Accept-Language", Language);

		StringBuilder sb = new StringBuilder();
		for (String sequence : sequences) {
			sb.append("input[]=").append(sequence).append("&");
		}
		sb.deleteCharAt(sb.length() - 1);

		String urlParameters = sb.toString();

		// Send post request
		con.setDoOutput(true);
		con.setDoInput(true);

		try {
			DataOutputStream wr = new DataOutputStream(con.getOutputStream());
			wr.writeBytes(urlParameters);
			wr.flush();
			wr.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			LOGGER.error("Error in matching peptides to LCA by " + pept2lca, e);
		}

		HashMap<String, UnipLCAResult> map = new HashMap<String, UnipLCAResult>();

		String line = null;
		try {
			BufferedReader in = new BufferedReader(new InputStreamReader(con.getInputStream()));
			line = in.readLine();
			in.close();
			LOGGER.debug("Received response from Unipept API, response length: " + (line != null ? line.length() : 0));
		} catch (IOException e) {
			// TODO Auto-generated catch block
			LOGGER.error("Error in matching peptides to LCA by " + pept2lca, e);
		}

		if (line == null || line.length() == 0) {
			LOGGER.warn("Empty response from Unipept API for " + sequences.length + " sequences");
			return map;
		}

		String[] cs = line.split("\\},\\{");
		LOGGER.debug("Parsing " + cs.length + " response entries from Unipept API");
		
		int successCount = 0;
		int errorCount = 0;
		
		for (int i = 0; i < cs.length; i++) {
			String[] csi = cs[i].split("[,:\"]+");
			if (csi.length == 9) {
				UnipLCAResult result = new UnipLCAResult(csi[2], Integer.parseInt(csi[4]), csi[6], csi[8]);
				map.put(result.getSequence(), result);
				successCount++;
			} else if (csi.length == 10) {
				UnipLCAResult result = new UnipLCAResult(csi[2], Integer.parseInt(csi[4]), csi[6], csi[8]);
				map.put(result.getSequence(), result);
				successCount++;
			} else {
				LOGGER.warn("Unexpected response format for entry " + i + ": expected 9 or 10 fields, got " + csi.length + ". Entry: " + cs[i]);
				errorCount++;
			}
		}

		LOGGER.info("pept2lca completed: " + successCount + " successful mappings, " + errorCount + " errors, " + 
				   "final map size: " + map.size() + " out of " + sequences.length + " input sequences");
		
		return map;
	}

	public static HashMap<String, UnipProInfo> protinfo(String[] pro) {

		URL obj = null;
		try {
			obj = new URL(protinfo);
		} catch (MalformedURLException e) {
			// TODO Auto-generated catch block
			LOGGER.error("Error in getting protein information by " + protinfo, e);
			return null;
		}
		HttpURLConnection con = null;
		try {
			con = (HttpURLConnection) obj.openConnection();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			LOGGER.error("Error in getting protein information by " + protinfo, e);
			return null;
		}

		// add request header
		try {
			con.setRequestMethod("POST");
		} catch (ProtocolException e) {
			// TODO Auto-generated catch block
			LOGGER.error("Error in getting protein information by " + protinfo, e);
		}
		con.setRequestProperty("accept", "*/*");
		con.setRequestProperty("connection", "Keep-Alive");
		con.setRequestProperty("User-Agent", USER_AGENT);
		con.setRequestProperty("Accept-Language", Language);

		StringBuilder sb = new StringBuilder();
		for (String protein : pro) {
			sb.append("input[]=").append(protein).append("&");
		}
		sb.deleteCharAt(sb.length() - 1);

		String urlParameters = sb.toString();

		// Send post request
		con.setDoOutput(true);
		con.setDoInput(true);

		try {
			DataOutputStream wr = new DataOutputStream(con.getOutputStream());
			wr.writeBytes(urlParameters);
			wr.flush();
			wr.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			LOGGER.error("Error in getting protein information by " + protinfo, e);
		}

		HashMap<String, UnipProInfo> map = new HashMap<String, UnipProInfo>();

		String line = null;
		try {
			BufferedReader in = new BufferedReader(new InputStreamReader(con.getInputStream()));
			line = in.readLine();
			in.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			LOGGER.error("Error in getting protein information by " + protinfo, e);
		}

		if (line == null || line.length() < 4) {
			return map;
		}

		try {
			// Remove the outer brackets and split by },{
			String jsonContent = line.substring(1, line.length() - 1);
			String[] entries = jsonContent.split("\\},\\{");

			for (int i = 0; i < entries.length; i++) {
				String entry = entries[i];
				
				// Parse JSON-like structure
				// Expected format: {"protein":"...","name":"...","taxon_id":...,"ec":[...],"go":[...],"ipr":[...]}
				String protein = extractJsonValue(entry, "protein");
				String name = extractJsonValue(entry, "name");
				String taxonIdStr = extractJsonValue(entry, "taxon_id");
				
				// Parse EC numbers
				String[] ec = extractArrayValues(entry, "ec");
				
				// Parse GO terms
				String[] go = extractArrayValues(entry, "go");
				
				if (!protein.isEmpty()) {
					int taxonId = 0;
					try {
						taxonId = Integer.parseInt(taxonIdStr);
					} catch (NumberFormatException e) {
						LOGGER.warn("Error parsing taxon ID '" + taxonIdStr + "': " + e.getMessage());
					}
					
					UnipProInfo result = new UnipProInfo(protein, name, taxonId, ec, go);
					map.put(protein, result);
				}
			}
			
		} catch (Exception e) {
			LOGGER.error("Error in parsing protinfo response: " + e.getMessage(), e);
		}

		return map;
	}

	/**
	 * Extract array values from JSON-like string (for EC, GO, IPR arrays)
	 */
	private static String[] extractArrayValues(String json, String key) {
		try {
			// Look for the key followed by an array
			String pattern = "\"" + key + "\"\\s*:\\s*\\[([^\\]]*)\\]";
			java.util.regex.Pattern p = java.util.regex.Pattern.compile(pattern);
			java.util.regex.Matcher m = p.matcher(json);
			
			if (m.find()) {
				String arrayContent = m.group(1).trim();
				if (arrayContent.isEmpty()) {
					return new String[0];
				}
				
				// Split by },{ to get individual objects
				String[] items = arrayContent.split("\\},\\{");
				ArrayList<String> values = new ArrayList<String>();
				
				for (String item : items) {
					// Extract the value from each object
					// For EC: {"ec_number": "1.1.1.1"} -> "1.1.1.1"
					// For GO: {"go_term": "GO:0032421"} -> "GO:0032421"
					// For IPR: {"code": "IPR026748"} -> "IPR026748"
					
					if (key.equals("ec")) {
						String ecNumber = extractJsonValue(item, "ec_number");
						if (!ecNumber.isEmpty()) {
							values.add(ecNumber);
						}
					} else if (key.equals("go")) {
						String goTerm = extractJsonValue(item, "go_term");
						if (!goTerm.isEmpty()) {
							values.add(goTerm);
						}
					} else if (key.equals("ipr")) {
						String iprCode = extractJsonValue(item, "code");
						if (!iprCode.isEmpty()) {
							values.add(iprCode);
						}
					}
				}
				
				return values.toArray(new String[values.size()]);
			}
			
			return new String[0];
		} catch (Exception e) {
			LOGGER.error("Error extracting array values for key '" + key + "': " + e.getMessage());
			return new String[0];
		}
	}

	/**
	 * Extract a value from a JSON-like string
	 */
	private static String extractJsonValue(String json, String key) {
		try {
			String pattern = "\"" + key + "\"\\s*:\\s*\"([^\"]+)\"";
			java.util.regex.Pattern p = java.util.regex.Pattern.compile(pattern);
			java.util.regex.Matcher m = p.matcher(json);
			
			if (m.find()) {
				return m.group(1);
			}
			
			// Try without quotes for numbers
			pattern = "\"" + key + "\"\\s*:\\s*(\\d+)";
			p = java.util.regex.Pattern.compile(pattern);
			m = p.matcher(json);
			
			if (m.find()) {
				return m.group(1);
			}
			
			return "";
		} catch (Exception e) {
			LOGGER.error("Error extracting JSON value for key '" + key + "': " + e.getMessage());
			return "";
		}
	}

	@SuppressWarnings("unused")
	private static void test(String in) {
		try (BufferedReader reader = new BufferedReader(new FileReader(in))) {
			String line = reader.readLine();
			ArrayList<String> list = new ArrayList<String>();

			while ((line = reader.readLine()) != null) {
				String[] cs = line.split("\t");
				list.add(cs[0]);
			}
			reader.close();

			String[] sequences = list.toArray(String[]::new);
			HashMap<String, UnipLCAResult> lcaResultMap = UnipeptRequester.pept2lca(sequences);
			System.out.println(sequences.length + "\t" + lcaResultMap.size());

		} catch (IOException e) {

		}
	}
}
