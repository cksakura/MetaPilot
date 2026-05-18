package uk.ac.quadram.metapilot.core.runner;

import uk.ac.quadram.metapilot.core.taxonomy.TaxonomyRanks;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;

/**
 * Static Genome.html writer shared by GUI and headless workflows.
 */
public final class GenomeHtmlWriter {

	private static final String LINK = "https://www.ebi.ac.uk/metagenomics/genomes/";

	private static final String HEAD = "<!DOCTYPE html>\n<html lang=\"en\">\n<head>\n"
			+ "<meta charset=\"UTF-8\">\n"
			+ "<meta name=\"viewport\" content=\"width=device-width, initial-scale=1.0\">\n"
			+ "<title>MetaPilot - Genome Analysis</title>\n"
			+ "<link rel=\"shortcut icon\" type=\"image/x-icon\" href=\"resources/favicon.ico\"/>\n"
			+ "<link href=\"https://cdnjs.cloudflare.com/ajax/libs/font-awesome/6.0.0/css/all.min.css\" rel=\"stylesheet\">\n"
			+ "<style>\n"
			+ "* { box-sizing: border-box; margin: 0; padding: 0; }\n"
			+ "body { font-family: 'Segoe UI', 'Roboto', 'Helvetica Neue', Arial, sans-serif; background: linear-gradient(135deg, #f8f9fa 0%, #e9ecef 100%); color: #333; line-height: 1.6; min-height: 100vh; }\n"
			+ ".container { max-width: 1400px; margin: 0 auto; padding: 20px; }\n"
			+ ".header { background: linear-gradient(135deg, #B5BD00 0%, #8ED1FC 100%); color: #000; padding: 30px; border-radius: 10px; margin-bottom: 30px; box-shadow: 0 4px 15px rgba(0,0,0,0.1); text-align: center; }\n"
			+ ".header h1 { font-size: 2.5em; margin-bottom: 10px; font-weight: 700; }\n"
			+ ".header p { font-size: 1.2em; opacity: 0.9; }\n"
			+ ".stats { display: flex; justify-content: space-around; margin: 20px 0; flex-wrap: wrap; }\n"
			+ ".stat-item { background: white; padding: 20px; border-radius: 8px; box-shadow: 0 2px 10px rgba(0,0,0,0.1); text-align: center; margin: 10px; min-width: 150px; }\n"
			+ ".stat-number { font-size: 2em; font-weight: bold; color: #B5BD00; }\n"
			+ ".stat-label { color: #666; margin-top: 5px; }\n"
			+ ".table-container { background: white; border-radius: 10px; box-shadow: 0 4px 15px rgba(0,0,0,0.1); overflow: hidden; margin-bottom: 30px; }\n"
			+ "table { width: 100%; border-collapse: collapse; font-size: 0.9em; }\n"
			+ "th { background: linear-gradient(135deg, #B5BD00 0%, #8ED1FC 100%); color: #000; padding: 15px 10px; text-align: left; font-weight: 600; position: sticky; top: 0; z-index: 10; }\n"
			+ "td { padding: 12px 10px; border-bottom: 1px solid #eee; vertical-align: top; }\n"
			+ "tr:hover { background-color: #f8f9fa; }\n"
			+ "tr:nth-child(even) { background-color: #f8f9fa; }\n"
			+ "tr:nth-child(even):hover { background-color: #e9ecef; }\n"
			+ "a { color: #B5BD00; text-decoration: none; font-weight: 500; transition: color 0.3s ease; }\n"
			+ "a:hover { color: #8ED1FC; text-decoration: underline; }\n"
			+ ".genome-link { font-family: 'Courier New', monospace; background: #f8f9fa; padding: 4px 8px; border-radius: 4px; border: 1px solid #dee2e6; }\n"
			+ ".taxonomy-cell { max-width: 200px; word-wrap: break-word; }\n"
			+ ".rank-badge { background: #8ED1FC; color: #000; padding: 2px 8px; border-radius: 12px; font-size: 0.8em; font-weight: 500; }\n"
			+ ".footer { text-align: center; padding: 20px; color: #666; background: white; border-radius: 10px; box-shadow: 0 2px 10px rgba(0,0,0,0.1); }\n"
			+ "@media (max-width: 768px) { .container { padding: 10px; } .header h1 { font-size: 2em; } .stats { flex-direction: column; } table { font-size: 0.8em; } th, td { padding: 8px 5px; } }\n"
			+ "</style>\n</head>\n<body>\n"
			+ "<div class=\"container\">\n"
			+ "<div class=\"header\">\n"
			+ "<h1><i class=\"fas fa-dna\"></i> MetaPilot Genome Analysis</h1>\n"
			+ "<p>Comprehensive taxonomic classification of identified genomes</p>\n"
			+ "</div>\n";

	private final PrintWriter writer;

	public GenomeHtmlWriter(File output) throws IOException {
		this.writer = new PrintWriter(output);
		this.writer.write(HEAD);
	}

	public void write(HashMap<String, String[]> genomeMap) {
		String[] genomes = genomeMap.keySet().toArray(new String[0]);
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

		writer.println("<div class=\"stats\">");
		writer.println("<div class=\"stat-item\"><div class=\"stat-number\">" + genomes.length
				+ "</div><div class=\"stat-label\">Total Genomes</div></div>");
		for (TaxonomyRanks rank : mainRanks) {
			String rankName = rank.getName();
			if(!rankName.equalsIgnoreCase("Kingdom")){
				writer.println("<div class=\"stat-item\"><div class=\"stat-number\">"
						+ rankTaxaCounts.get(rankName).size() + "</div><div class=\"stat-label\">"
						+ formatRankLabel(rankName) + "</div></div>");
			}
		}
		writer.println("</div>");
		writer.println("<div class=\"table-container\"><table><tr>");

		ArrayList<String> list = new ArrayList<>();
		list.add("Genome ID");
		list.add("Taxonomy");
		list.add("Rank");
		for (String title : list) {
			writer.println("<th>" + title + "</th>");
		}
		for (TaxonomyRanks rank : mainRanks) {
			writer.println("<th>" + rank.getName() + "</th>");
		}
		writer.println("</tr>");

		for (String genome : genomes) {
			String[] taxa = genomeMap.get(genome);
			writer.println("<tr>");
			writer.println("<td><a href=\"" + LINK + genome + "\" class=\"genome-link\" target=\"_blank\">"
					+ genome + "</a></td>");
			String taxonomyName = taxa.length > 1 ? taxa[1] : "";
			writer.println("<td class=\"taxonomy-cell\">" + taxonomyName + "</td>");
			String rank = taxa.length > 2 ? taxa[2] : "";
			writer.println("<td><span class=\"rank-badge\">" + rank + "</span></td>");
			for (int i = 3; i < mainRanks.length + 3 && i < taxa.length; i++) {
				String taxonValue = taxa[i] != null ? taxa[i] : "";
				writer.println("<td class=\"taxonomy-cell\">" + taxonValue + "</td>");
			}
			writer.println("</tr>");
		}
		writer.println("</table></div>");
	}

	public void close() {
		writer.println("<div class=\"footer\">");
		writer.println("<p><i class=\"fas fa-info-circle\"></i> Generated by MetaPilot Genome Analysis</p>");
		writer.println("<p>Click on genome IDs to view detailed information in the EBI Metagenomics database</p>");
		writer.println("</div>");
		writer.println("</div>");
		writer.println("</body>");
		writer.println("</html>");
		writer.close();
	}

	private static String formatRankLabel(String rankName) {
		if (rankName == null || rankName.isEmpty()) {
			return "";
		}
		return Character.toUpperCase(rankName.charAt(0)) + rankName.substring(1);
	}
}
