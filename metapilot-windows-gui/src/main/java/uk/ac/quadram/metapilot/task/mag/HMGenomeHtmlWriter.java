package uk.ac.quadram.metapilot.task.mag;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;

import uk.ac.quadram.metapilot.core.taxonomy.TaxonomyRanks;

public class HMGenomeHtmlWriter {

	private String link = "https://www.ebi.ac.uk/metagenomics/genomes/";

	private static final String head = "<!DOCTYPE html>\n<html lang=\"en\">\n<head>\n"
			+ "<meta charset=\"UTF-8\">\n"
			+ "<meta name=\"viewport\" content=\"width=device-width, initial-scale=1.0\">\n"
			+ "<title>MetaPilot - Genome Analysis</title>\n"
			+ "<link rel=\"shortcut icon\" type=\"image/x-icon\" href=\"resources/favicon.ico\"/>\n"
			+ "<link href=\"https://cdnjs.cloudflare.com/ajax/libs/font-awesome/6.0.0/css/all.min.css\" rel=\"stylesheet\">\n"
			+ "<style>\n"
			+ "/* MetaPilot Genome Report Styling */\n"
			+ "* {\n"
			+ "    box-sizing: border-box;\n"
			+ "    margin: 0;\n"
			+ "    padding: 0;\n"
			+ "}\n\n"
			+ "body {\n"
			+ "    font-family: 'Segoe UI', 'Roboto', 'Helvetica Neue', Arial, sans-serif;\n"
			+ "    background: linear-gradient(135deg, #f8f9fa 0%, #e9ecef 100%);\n"
			+ "    color: #333;\n"
			+ "    line-height: 1.6;\n"
			+ "    min-height: 100vh;\n"
			+ "}\n\n"
			+ ".container {\n"
			+ "    max-width: 1400px;\n"
			+ "    margin: 0 auto;\n"
			+ "    padding: 20px;\n"
			+ "}\n\n"
			+ ".header {\n"
			+ "    background: linear-gradient(135deg, #B5BD00 0%, #8ED1FC 100%);\n"
			+ "    color: #000;\n"
			+ "    padding: 30px;\n"
			+ "    border-radius: 10px;\n"
			+ "    margin-bottom: 30px;\n"
			+ "    box-shadow: 0 4px 15px rgba(0,0,0,0.1);\n"
			+ "    text-align: center;\n"
			+ "}\n\n"
			+ ".header h1 {\n"
			+ "    font-size: 2.5em;\n"
			+ "    margin-bottom: 10px;\n"
			+ "    font-weight: 700;\n"
			+ "}\n\n"
			+ ".header p {\n"
			+ "    font-size: 1.2em;\n"
			+ "    opacity: 0.9;\n"
			+ "}\n\n"
			+ ".stats {\n"
			+ "    display: flex;\n"
			+ "    justify-content: space-around;\n"
			+ "    margin: 20px 0;\n"
			+ "    flex-wrap: wrap;\n"
			+ "}\n\n"
			+ ".stat-item {\n"
			+ "    background: white;\n"
			+ "    padding: 20px;\n"
			+ "    border-radius: 8px;\n"
			+ "    box-shadow: 0 2px 10px rgba(0,0,0,0.1);\n"
			+ "    text-align: center;\n"
			+ "    margin: 10px;\n"
			+ "    min-width: 150px;\n"
			+ "}\n\n"
			+ ".stat-number {\n"
			+ "    font-size: 2em;\n"
			+ "    font-weight: bold;\n"
			+ "    color: #B5BD00;\n"
			+ "}\n\n"
			+ ".stat-label {\n"
			+ "    color: #666;\n"
			+ "    margin-top: 5px;\n"
			+ "}\n\n"
			+ ".table-container {\n"
			+ "    background: white;\n"
			+ "    border-radius: 10px;\n"
			+ "    box-shadow: 0 4px 15px rgba(0,0,0,0.1);\n"
			+ "    overflow: hidden;\n"
			+ "    margin-bottom: 30px;\n"
			+ "}\n\n"
			+ "table {\n"
			+ "    width: 100%;\n"
			+ "    border-collapse: collapse;\n"
			+ "    font-size: 0.9em;\n"
			+ "}\n\n"
			+ "th {\n"
			+ "    background: linear-gradient(135deg, #B5BD00 0%, #8ED1FC 100%);\n"
			+ "    color: #000;\n"
			+ "    padding: 15px 10px;\n"
			+ "    text-align: left;\n"
			+ "    font-weight: 600;\n"
			+ "    position: sticky;\n"
			+ "    top: 0;\n"
			+ "    z-index: 10;\n"
			+ "}\n\n"
			+ "td {\n"
			+ "    padding: 12px 10px;\n"
			+ "    border-bottom: 1px solid #eee;\n"
			+ "    vertical-align: top;\n"
			+ "}\n\n"
			+ "tr:hover {\n"
			+ "    background-color: #f8f9fa;\n"
			+ "}\n\n"
			+ "tr:nth-child(even) {\n"
			+ "    background-color: #f8f9fa;\n"
			+ "}\n\n"
			+ "tr:nth-child(even):hover {\n"
			+ "    background-color: #e9ecef;\n"
			+ "}\n\n"
			+ "a {\n"
			+ "    color: #B5BD00;\n"
			+ "    text-decoration: none;\n"
			+ "    font-weight: 500;\n"
			+ "    transition: color 0.3s ease;\n"
			+ "}\n\n"
			+ "a:hover {\n"
			+ "    color: #8ED1FC;\n"
			+ "    text-decoration: underline;\n"
			+ "}\n\n"
			+ ".genome-link {\n"
			+ "    font-family: 'Courier New', monospace;\n"
			+ "    background: #f8f9fa;\n"
			+ "    padding: 4px 8px;\n"
			+ "    border-radius: 4px;\n"
			+ "    border: 1px solid #dee2e6;\n"
			+ "}\n\n"
			+ ".taxonomy-cell {\n"
			+ "    max-width: 200px;\n"
			+ "    word-wrap: break-word;\n"
			+ "}\n\n"
			+ ".rank-badge {\n"
			+ "    background: #8ED1FC;\n"
			+ "    color: #000;\n"
			+ "    padding: 2px 8px;\n"
			+ "    border-radius: 12px;\n"
			+ "    font-size: 0.8em;\n"
			+ "    font-weight: 500;\n"
			+ "}\n\n"
			+ ".footer {\n"
			+ "    text-align: center;\n"
			+ "    padding: 20px;\n"
			+ "    color: #666;\n"
			+ "    background: white;\n"
			+ "    border-radius: 10px;\n"
			+ "    box-shadow: 0 2px 10px rgba(0,0,0,0.1);\n"
			+ "}\n\n"
			+ "@media (max-width: 768px) {\n"
			+ "    .container {\n"
			+ "        padding: 10px;\n"
			+ "    }\n"
			+ "    .header h1 {\n"
			+ "        font-size: 2em;\n"
			+ "    }\n"
			+ "    .stats {\n"
			+ "        flex-direction: column;\n"
			+ "    }\n"
			+ "    table {\n"
			+ "        font-size: 0.8em;\n"
			+ "    }\n"
			+ "    th, td {\n"
			+ "        padding: 8px 5px;\n"
			+ "    }\n"
			+ "}\n"
			+ "</style>\n</head>\n<body>\n"
			+ "<div class=\"container\">\n"
			+ "<div class=\"header\">\n"
			+ "<h1><i class=\"fas fa-dna\"></i> MetaPilot Genome Analysis</h1>\n"
			+ "<p>Comprehensive taxonomic classification of identified genomes</p>\n"
			+ "</div>\n";

	private PrintWriter writer;

	public HMGenomeHtmlWriter(String output) throws IOException {
		this(new File(output));
	}

	public HMGenomeHtmlWriter(File output) throws IOException {
		this.writer = new PrintWriter(output);
		this.writer.write(head);
	}

	public void write(HashMap<String, String[]> genomeMap) {
		// Generate statistics
		String[] genomes = genomeMap.keySet().toArray(new String[genomeMap.size()]);
		Arrays.sort(genomes);
		
		TaxonomyRanks[] mainRanks = TaxonomyRanks.getMainRanks();
		HashMap<String, HashSet<String>> rankTaxaCounts = new HashMap<>();
		for (TaxonomyRanks rank : mainRanks) {
			rankTaxaCounts.put(rank.getName(), new HashSet<>());
		}
		
		for (String genome : genomes) {
			String[] taxa = genomeMap.get(genome);
			for (int i = 0; i < mainRanks.length; i++) {
				int taxaIndex = i + 3;
				if (taxa.length > taxaIndex) {
					String taxon = taxa[taxaIndex];
					if (taxon != null && !taxon.trim().isEmpty()) {
						rankTaxaCounts.get(mainRanks[i].getName()).add(taxon.trim());
					}
				}
			}
		}
		
		// Add statistics section
		this.writer.println("<div class=\"stats\">");
		this.writer.println("<div class=\"stat-item\">");
		this.writer.println("<div class=\"stat-number\">" + genomes.length + "</div>");
		this.writer.println("<div class=\"stat-label\">Total Genomes</div>");
		this.writer.println("</div>");
		for (TaxonomyRanks rank : mainRanks) {
			String rankName = rank.getName();
			this.writer.println("<div class=\"stat-item\">");
			this.writer.println("<div class=\"stat-number\">" + rankTaxaCounts.get(rankName).size() + "</div>");
			this.writer.println("<div class=\"stat-label\">" + formatRankLabel(rankName) + "</div>");
			this.writer.println("</div>");
		}
		this.writer.println("</div>");

		// Add table container
		this.writer.println("<div class=\"table-container\">");
		this.writer.println("<table>");
		this.writer.println("<tr>");

		ArrayList<String> list = new ArrayList<String>();
		list.add("Genome ID");
		list.add("Taxonomy");
		list.add("Rank");

		for (String title : list) {
			this.writer.println("<th>" + title + "</th>");
		}
		for (TaxonomyRanks rank : mainRanks) {
			this.writer.println("<th>" + rank.getName() + "</th>");
		}
		this.writer.println("</tr>");

		for (String genome : genomes) {
			String[] taxa = genomeMap.get(genome);
			this.writer.println("<tr>");
			
			// Genome ID with enhanced styling
			this.writer.println("<td><a href=\"" + link + genome + "\" class=\"genome-link\" target=\"_blank\">" + genome + "</a></td>");
			
			// Taxonomy name
			String taxonomyName = taxa.length > 1 ? taxa[1] : "";
			this.writer.println("<td class=\"taxonomy-cell\">" + taxonomyName + "</td>");
			
			// Rank with badge styling
			String rank = taxa.length > 2 ? taxa[2] : "";
			this.writer.println("<td><span class=\"rank-badge\">" + rank + "</span></td>");
			
			// Taxonomic ranks
			for (int i = 3; i < mainRanks.length + 3 && i < taxa.length; i++) {
				String taxonValue = taxa[i] != null ? taxa[i] : "";
				this.writer.println("<td class=\"taxonomy-cell\">" + taxonValue + "</td>");
			}
			this.writer.println("</tr>");
		}
		this.writer.println("</table>");
		this.writer.println("</div>");
	}

	public void close() {
		this.writer.println("<div class=\"footer\">");
		this.writer.println("<p><i class=\"fas fa-info-circle\"></i> Generated by MetaPilot Genome Analysis</p>");
		this.writer.println("<p>Click on genome IDs to view detailed information in the EBI Metagenomics database</p>");
		this.writer.println("</div>");
		this.writer.println("</div>"); // Close container
		this.writer.println("</body>");
		this.writer.println("</html>");
		this.writer.close();
	}

	private static String formatRankLabel(String rankName) {
		if (rankName == null || rankName.isEmpty()) {
			return "";
		}
		return Character.toUpperCase(rankName.charAt(0)) + rankName.substring(1);
	}
}
