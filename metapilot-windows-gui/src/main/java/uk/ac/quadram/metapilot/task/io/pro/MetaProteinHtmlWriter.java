/**
 * 
 */
package uk.ac.quadram.metapilot.task.io.pro;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.time.LocalDateTime;
import java.time.format.DateTimeFormatter;
import java.util.ArrayList;
import java.util.HashMap;

/**
 * @author Kai Cheng
 *
 */
public class MetaProteinHtmlWriter {

	private String Gene_Ontology_link = "https://www.ebi.ac.uk/QuickGO/term/";
	private String EC_link = "https://www.genome.jp/dbget-bin/www_bget?ec:";
	private String KEGG_link = "https://www.genome.jp/dbget-bin/www_bget?";
	private String BRITE_link = "https://www.genome.jp/kegg/annotation/";
	private String KEGG_TC_link = "http://www.tcdb.org/search/result.php?tc=";
	private String CAZy_link = "http://www.cazy.org/";
	private String BiGG_Reaction_link = "http://bigg.ucsd.edu/models/";

	// Helper methods for generating links (replacing the removed Finder classes)
	
	/**
	 * Generate COG database link
	 */
	private String getCOGLink(String cog) {
		if (cog == null || cog.trim().isEmpty()) {
			return "";
		}
		return "https://www.ncbi.nlm.nih.gov/research/cog/" + cog;
	}
	
	/**
	 * Generate KEGG link for KO or pathway entries
	 */
	private String getKEGGLink(String kegg) {
		if (kegg == null || kegg.trim().isEmpty()) {
			return "";
		}
		// Handle different KEGG entry types
		if (kegg.startsWith("ko:")) {
			return "https://www.genome.jp/dbget-bin/www_bget?" + kegg;
		} else if (kegg.startsWith("map")) {
			return "https://www.genome.jp/dbget-bin/www_bget?path:" + kegg;
		} else if (kegg.startsWith("M")) {
			return "https://www.genome.jp/module/" + kegg;
		} else if (kegg.startsWith("R")) {
			return "https://www.genome.jp/dbget-bin/www_bget?rn:" + kegg;
		} else {
			return "https://www.genome.jp/dbget-bin/www_bget?" + kegg;
		}
	}
	
	/**
	 * Generate QuickGO link for GO Biological Process terms
	 */
	private String getGOBPQuickGoLink(String go) {
		if (go == null || go.trim().isEmpty()) {
			return "";
		}
		return "https://www.ebi.ac.uk/QuickGO/term/" + go;
	}
	
	/**
	 * Generate AmiGO link for GO Biological Process terms
	 */
	private String getGOBPAmiGoLink(String go) {
		if (go == null || go.trim().isEmpty()) {
			return "";
		}
		return "http://amigo.geneontology.org/amigo/term/" + go;
	}
	
	/**
	 * Generate QuickGO link for GO Cellular Component terms
	 */
	private String getGOCCQuickGoLink(String go) {
		if (go == null || go.trim().isEmpty()) {
			return "";
		}
		return "https://www.ebi.ac.uk/QuickGO/term/" + go;
	}
	
	/**
	 * Generate AmiGO link for GO Cellular Component terms
	 */
	private String getGOCCAmiGoLink(String go) {
		if (go == null || go.trim().isEmpty()) {
			return "";
		}
		return "http://amigo.geneontology.org/amigo/term/" + go;
	}
	
	/**
	 * Generate QuickGO link for GO Molecular Function terms
	 */
	private String getGOMFQuickGoLink(String go) {
		if (go == null || go.trim().isEmpty()) {
			return "";
		}
		return "https://www.ebi.ac.uk/QuickGO/term/" + go;
	}
	
	/**
	 * Generate AmiGO link for GO Molecular Function terms
	 */
	private String getGOMFAmiGoLink(String go) {
		if (go == null || go.trim().isEmpty()) {
			return "";
		}
		return "http://amigo.geneontology.org/amigo/term/" + go;
	}

	private static final String head = "<!DOCTYPE html>\n<html lang=\"en\">\n<head>\n"
			+ "<meta charset=\"UTF-8\">\n"
			+ "<meta name=\"viewport\" content=\"width=device-width, initial-scale=1.0\">\n"
			+ "<title>MetaPilot - Protein Functional Analysis</title>\n"
			+ "<link rel=\"shortcut icon\" type=\"image/x-icon\" href=\"resources/favicon.ico\"/>\n"
			+ "<link href=\"https://cdnjs.cloudflare.com/ajax/libs/font-awesome/6.0.0/css/all.min.css\" rel=\"stylesheet\">\n"
			+ "<style>\n"
			+ "/* MetaPilot Protein Report Styling */\n"
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
			+ "    box-shadow: 0 8px 25px rgba(181, 189, 0, 0.2);\n"
			+ "    text-align: center;\n"
			+ "}\n\n"
			+ ".header h1 {\n"
			+ "    font-size: 2.5rem;\n"
			+ "    font-weight: 700;\n"
			+ "    margin-bottom: 10px;\n"
			+ "    text-shadow: 2px 2px 4px rgba(0,0,0,0.1);\n"
			+ "}\n\n"
			+ ".header p {\n"
			+ "    font-size: 1.1rem;\n"
			+ "    opacity: 0.8;\n"
			+ "}\n\n"
			+ ".content-wrapper {\n"
			+ "    background: white;\n"
			+ "    border-radius: 10px;\n"
			+ "    box-shadow: 0 5px 20px rgba(0,0,0,0.1);\n"
			+ "    overflow: hidden;\n"
			+ "}\n\n"
			+ ".table-header {\n"
			+ "    background: linear-gradient(135deg, #2c3e50 0%, #34495e 100%);\n"
			+ "    color: white;\n"
			+ "    padding: 20px;\n"
			+ "    text-align: center;\n"
			+ "}\n\n"
			+ ".table-container {\n"
			+ "    overflow-x: auto;\n"
			+ "    padding: 0;\n"
			+ "}\n\n"
			+ "table {\n"
			+ "    width: 100%;\n"
			+ "    border-collapse: collapse;\n"
			+ "    font-size: 14px;\n"
			+ "    background: white;\n"
			+ "}\n\n"
			+ "th {\n"
			+ "    background: linear-gradient(135deg, #B5BD00 0%, #A4AC00 100%);\n"
			+ "    color: #000;\n"
			+ "    font-weight: 600;\n"
			+ "    padding: 15px 12px;\n"
			+ "    text-align: left;\n"
			+ "    position: sticky;\n"
			+ "    top: 0;\n"
			+ "    z-index: 10;\n"
			+ "    border-bottom: 2px solid #8ED1FC;\n"
			+ "    font-size: 13px;\n"
			+ "    text-transform: uppercase;\n"
			+ "    letter-spacing: 0.5px;\n"
			+ "}\n\n"
			+ "td {\n"
			+ "    padding: 12px;\n"
			+ "    border-bottom: 1px solid #e9ecef;\n"
			+ "    vertical-align: top;\n"
			+ "    transition: background-color 0.2s ease;\n"
			+ "}\n\n"
			+ "tr:nth-child(even) {\n"
			+ "    background-color: #f8f9fa;\n"
			+ "}\n\n"
			+ "tr:hover {\n"
			+ "    background-color: rgba(142, 209, 252, 0.1);\n"
			+ "    transform: translateX(2px);\n"
			+ "    transition: all 0.2s ease;\n"
			+ "}\n\n"
			+ "a {\n"
			+ "    color: #B5BD00;\n"
			+ "    text-decoration: none;\n"
			+ "    font-weight: 500;\n"
			+ "    transition: all 0.2s ease;\n"
			+ "    display: inline-block;\n"
			+ "    padding: 2px 6px;\n"
			+ "    border-radius: 4px;\n"
			+ "}\n\n"
			+ "a:hover {\n"
			+ "    color: #000;\n"
			+ "    background-color: #8ED1FC;\n"
			+ "    transform: translateY(-1px);\n"
			+ "    box-shadow: 0 2px 8px rgba(142, 209, 252, 0.3);\n"
			+ "}\n\n"
			+ "ul {\n"
			+ "    list-style: none;\n"
			+ "    padding: 0;\n"
			+ "    margin: 0;\n"
			+ "}\n\n"
			+ "li {\n"
			+ "    margin-bottom: 8px;\n"
			+ "    padding: 6px 10px;\n"
			+ "    background: rgba(181, 189, 0, 0.05);\n"
			+ "    border-left: 3px solid #B5BD00;\n"
			+ "    border-radius: 4px;\n"
			+ "    transition: all 0.2s ease;\n"
			+ "}\n\n"
			+ "li:hover {\n"
			+ "    background: rgba(142, 209, 252, 0.1);\n"
			+ "    border-left-color: #8ED1FC;\n"
			+ "    transform: translateX(3px);\n"
			+ "}\n\n"
			+ ".empty-cell {\n"
			+ "    color: #6c757d;\n"
			+ "    font-style: italic;\n"
			+ "    text-align: center;\n"
			+ "}\n\n"
			+ ".group-id {\n"
			+ "    font-weight: bold;\n"
			+ "    font-size: 1.1rem;\n"
			+ "    color: #2c3e50;\n"
			+ "    background: linear-gradient(135deg, #B5BD00, #8ED1FC);\n"
			+ "    -webkit-background-clip: text;\n"
			+ "    -webkit-text-fill-color: transparent;\n"
			+ "    background-clip: text;\n"
			+ "}\n\n"
			+ ".protein-names {\n"
			+ "    max-width: 200px;\n"
			+ "    word-wrap: break-word;\n"
			+ "}\n\n"
			+ ".footer {\n"
			+ "    text-align: center;\n"
			+ "    margin-top: 40px;\n"
			+ "    padding: 20px;\n"
			+ "    color: #6c757d;\n"
			+ "    background: white;\n"
			+ "    border-radius: 10px;\n"
			+ "    box-shadow: 0 2px 10px rgba(0,0,0,0.05);\n"
			+ "}\n\n"
			+ ".footer .logo {\n"
			+ "    color: #B5BD00;\n"
			+ "    font-weight: bold;\n"
			+ "    font-size: 1.2rem;\n"
			+ "}\n\n"
			+ "@media (max-width: 768px) {\n"
			+ "    .container { padding: 10px; }\n"
			+ "    .header h1 { font-size: 2rem; }\n"
			+ "    th, td { padding: 8px 6px; font-size: 12px; }\n"
			+ "    .protein-names { max-width: 150px; }\n"
			+ "}\n\n"
			+ ".icon {\n"
			+ "    margin-right: 8px;\n"
			+ "    color: #B5BD00;\n"
			+ "}\n"
			+ "</style>\n"
			+ "</head>\n"
			+ "<body>\n"
			+ "<div class=\"container\">\n"
			+ "    <div class=\"header\">\n"
			+ "        <h1><i class=\"fas fa-dna icon\"></i>MetaPilot Protein Functional Analysis</h1>\n"
			+ "        <p><i class=\"fas fa-chart-bar icon\"></i>Comprehensive protein annotation and functional classification</p>\n"
			+ "    </div>\n"
			+ "    <div class=\"content-wrapper\">\n"
			+ "        <div class=\"table-header\">\n"
			+ "            <h2><i class=\"fas fa-table icon\"></i>Protein Groups and Functional Annotations</h2>\n"
			+ "        </div>\n"
			+ "        <div class=\"table-container\">\n";

	private PrintWriter writer;

	public MetaProteinHtmlWriter(String output) throws IOException {
		this(new File(output));
	}

	public MetaProteinHtmlWriter(File output) throws IOException {
		this.writer = new PrintWriter(output);
		this.writer.write(head);
	}

	public void write(MetaProteinAnnoEggNog[] proteins, HashMap<String, String[]> cogMap,
			HashMap<String, String[]> nogMap) {
		this.writer.println("<table>");
		this.writer.println("<tr>");

		ArrayList<String> list = new ArrayList<String>();
		list.add("Group ID");
		list.add("Protein names");

		for (int i = 0; i < MetaProteinAnnoEggNog.funcNames.length; i++) {
			list.add(MetaProteinAnnoEggNog.funcNames[i]);
		}

		for (String title : list) {
			this.writer.println("<th>" + title + "</th>");
		}
		this.writer.println("</tr>");

		MetaProteinAnnoEggNog delegrate = null;
		ArrayList<String> proNameList = new ArrayList<String>();

		for (int i = 0; i < proteins.length; i++) {

			MetaProtein pro = proteins[i].getPro();

			if (delegrate == null) {
				delegrate = proteins[i];
				proNameList.add(pro.getName());
			} else {

				if (delegrate.getPro().getGroupId() != -1
						&& delegrate.getPro().getGroupId() == proteins[i].getPro().getGroupId()) {

					proNameList.add(proteins[i].getPro().getName());

				} else {

					this.writer.println("<tr>");
					this.writer.println("<td class=\"group-id\">" + (delegrate.getPro().getGroupId()) + "</td>");

					StringBuilder namesb = new StringBuilder("<td class=\"protein-names\">");
					for (String name : proNameList) {
						namesb.append(name).append(";").append("<br>");
					}

					namesb.append("</td>");
					this.writer.println(namesb);

					String Gene_Ontology = delegrate.getGene_Ontology();
					if (Gene_Ontology == null) {
						this.writer.println("<td class=\"empty-cell\">-</td>");
					} else {
						String[] gos = Gene_Ontology.split(",");
						if (gos != null && gos.length > 0) {
							StringBuilder sb = new StringBuilder("<td><ul>");
							for (String go : gos) {
								sb.append("<li><a href=\"").append(this.Gene_Ontology_link + go).append("\">")
										.append(go).append("</a></li>");

							}
							sb.append("</ul></td>");
							this.writer.println(sb);

						} else {
							this.writer.println("<td class=\"empty-cell\">-</td>");
						}
					}

					String EC = delegrate.getEC();
					if (EC == null) {
						this.writer.println("<td class=\"empty-cell\">-</td>");
					} else {
						String[] ecs = EC.split(",");
						if (ecs != null && ecs.length > 0) {
							StringBuilder sb = new StringBuilder("<td><ul>");
							for (String name : ecs) {
								sb.append("<li><a href=\"").append(this.EC_link + name).append("\">").append(name)
										.append("</a></li>");

							}
							sb.append("</ul></td>");
							this.writer.println(sb);

						} else {
							this.writer.println("<td class=\"empty-cell\">-</td>");
						}
					}

					String KEGG_ko = delegrate.getKEGG_ko();
					if (KEGG_ko == null) {
						this.writer.println("<td class=\"empty-cell\">-</td>");
					} else {
						String[] kos = KEGG_ko.split(",");
						if (kos != null && kos.length > 0) {
							StringBuilder sb = new StringBuilder("<td><ul>");
							for (String name : kos) {
								sb.append("<li><a href=\"").append(this.KEGG_link + name).append("\">").append(name)
										.append("</a></li>");

							}
							sb.append("</ul></td>");
							this.writer.println(sb);

						} else {
							this.writer.println("<td class=\"empty-cell\">-</td>");
						}
					}

					String KEGG_Pathway = delegrate.getKEGG_Pathway();
					if (KEGG_Pathway == null) {
						this.writer.println("<td class=\"empty-cell\">-</td>");
					} else {
						String[] pathways = KEGG_Pathway.split(",");
						if (pathways != null && pathways.length > 0) {
							StringBuilder sb = new StringBuilder("<td><ul>");
							for (String name : pathways) {
								sb.append("<li><a href=\"").append(this.KEGG_link + name).append("\">").append(name)
										.append("</a></li>");

							}
							sb.append("</ul></td>");
							this.writer.println(sb);

						} else {
							this.writer.println("<td class=\"empty-cell\">-</td>");
						}
					}

					String KEGG_Module = delegrate.getKEGG_Module();
					if (KEGG_Module == null) {
						this.writer.println("<td class=\"empty-cell\">-</td>");
					} else {
						String[] modules = KEGG_Module.split(",");
						if (modules != null && modules.length > 0) {
							StringBuilder sb = new StringBuilder("<td><ul>");
							for (String name : modules) {
								sb.append("<li><a href=\"").append(this.KEGG_link + name).append("\">").append(name)
										.append("</a></li>");

							}
							sb.append("</ul></td>");
							this.writer.println(sb);

						} else {
							this.writer.println("<td class=\"empty-cell\">-</td>");
						}
					}

					String KEGG_Reaction = delegrate.getKEGG_Reaction();
					if (KEGG_Reaction == null) {
						this.writer.println("<td class=\"empty-cell\">-</td>");
					} else {
						String[] reactions = KEGG_Reaction.split(",");
						if (reactions != null && reactions.length > 0) {
							StringBuilder sb = new StringBuilder("<td><ul>");
							for (String name : reactions) {
								String webname = name.indexOf(".") > 0 ? name.substring(0, name.indexOf(".")) : name;
								sb.append("<li><a href=\"").append(this.KEGG_link + webname).append("\">").append(name)
										.append("</a></li>");

							}
							sb.append("</ul></td>");
							this.writer.println(sb);

						} else {
							this.writer.println("<td class=\"empty-cell\">-</td>");
						}
					}

					String KEGG_rclass = delegrate.getKEGG_rclass();
					if (KEGG_rclass == null) {
						this.writer.println("<td class=\"empty-cell\">-</td>");
					} else {
						String[] rclass = KEGG_rclass.split(",");
						if (rclass != null && rclass.length > 0) {
							StringBuilder sb = new StringBuilder("<td><ul>");
							for (String name : rclass) {
								sb.append("<li><a href=\"").append(this.KEGG_link + name).append("\">").append(name)
										.append("</a></li>");

							}
							sb.append("</ul></td>");
							this.writer.println(sb);

						} else {
							this.writer.println("<td class=\"empty-cell\">-</td>");
						}
					}

					String BRITE = delegrate.getBRITE();
					if (BRITE == null) {
						this.writer.println("<td class=\"empty-cell\">-</td>");
					} else {
						String[] brites = BRITE.split(",");
						if (brites != null && brites.length > 0) {
							StringBuilder sb = new StringBuilder("<td><ul>");
							for (String name : brites) {
								if (name.startsWith("br")) {
									sb.append("<li><a href=\"").append(this.BRITE_link + name + ".html").append("\">")
											.append(name).append("</a></li>");
								} else {
									sb.append("<li>").append(name).append("</li>");
								}

							}
							sb.append("</ul></td>");
							this.writer.println(sb);

						} else {
							this.writer.println("<td class=\"empty-cell\">-</td>");
						}
					}

					String KEGG_TC = delegrate.getKEGG_TC();
					if (KEGG_TC == null) {
						this.writer.println("<td class=\"empty-cell\">-</td>");
					} else {
						String[] tcs = KEGG_TC.split(",");
						if (tcs != null && tcs.length > 0) {
							StringBuilder sb = new StringBuilder("<td><ul>");
							for (String name : tcs) {
								sb.append("<li><a href=\"").append(this.KEGG_TC_link + name).append("\">").append(name)
										.append("</a></li>");

							}
							sb.append("</ul></td>");
							this.writer.println(sb);

						} else {
							this.writer.println("<td class=\"empty-cell\">-</td>");
						}
					}

					String CAZy = delegrate.getCAZy();
					if (CAZy == null) {
						this.writer.println("<td class=\"empty-cell\">-</td>");
					} else {
						String[] cazys = CAZy.split(",");
						if (cazys != null && cazys.length > 0) {
							StringBuilder sb = new StringBuilder("<td><ul>");
							for (String name : cazys) {
								sb.append("<li><a href=\"").append(this.CAZy_link + name + ".html").append("\">")
										.append(name).append("</a></li>");

							}
							sb.append("</ul></td>");
							this.writer.println(sb);

						} else {
							this.writer.println("<td class=\"empty-cell\">-</td>");
						}
					}

					String BiGG_Reaction = delegrate.getBiGG_Reaction();
					if (BiGG_Reaction == null) {
						this.writer.println("<td class=\"empty-cell\">-</td>");
					} else {
						String[] biggs = BiGG_Reaction.split(",");
						if (biggs != null && biggs.length > 0) {
							StringBuilder sb = new StringBuilder("<td><ul>");
							for (String name : biggs) {
								if (name.indexOf(".") > 0) {
									sb.append("<li><a href=\"")
											.append(this.BiGG_Reaction_link + name.substring(0, name.indexOf(".")))
											.append("\">").append(name).append("</a></li>");
								} else {
									sb.append("<li><a href=\"").append(this.BiGG_Reaction_link + name).append("\">")
											.append(name).append("</a></li>");
								}
							}
							sb.append("</ul></td>");
							this.writer.println(sb);

						} else {
							this.writer.println("<td class=\"empty-cell\">-</td>");
						}
					}

					String cog = delegrate.getCog();
					if (cog == null) {
						this.writer.println("<td class=\"empty-cell\">-</td>");
						this.writer.println("<td class=\"empty-cell\">-</td>");
						this.writer.println("<td class=\"empty-cell\">-</td>");
					} else {
						this.writer.println("<td><a href=\"" + getCOGLink(cog) + "\">" + cog + "</a></td>");
						if (cogMap.containsKey(cog)) {
							this.writer.println("<td>" + cogMap.get(cog)[0] + "</td>");
							this.writer.println("<td>" + cogMap.get(cog)[1] + "</td>");
						} else {
							this.writer.println("<td class=\"empty-cell\">-</td>");
							this.writer.println("<td class=\"empty-cell\">-</td>");
						}
					}

					String nog = delegrate.getNog();
					if (nog == null) {
						this.writer.println("<td class=\"empty-cell\">-</td>");
						this.writer.println("<td class=\"empty-cell\">-</td>");
						this.writer.println("<td class=\"empty-cell\">-</td>");
					} else {
						this.writer.println("<td>" + nog + "</td>");
						if (nogMap.containsKey(nog)) {
							this.writer.println("<td>" + nogMap.get(nog)[0] + "</td>");
							this.writer.println("<td>" + nogMap.get(nog)[1] + "</td>");
						} else {
							this.writer.println("<td class=\"empty-cell\">-</td>");
							this.writer.println("<td class=\"empty-cell\">-</td>");
						}
					}

					this.writer.println("</tr>");

					delegrate = proteins[i];
					proNameList = new ArrayList<String>();
					proNameList.add(proteins[i].getPro().getName());
				}
			}
		}
	}

	public void write(MetaProteinAnno1[] proteins, HashMap<String, HashMap<String, String>> functionMap) {

		this.writer.println("<table>");
		this.writer.println("<tr>");

		ArrayList<String> list = new ArrayList<String>();
		list.add("Group ID");
		list.add("Protein names");

		HashMap<String, String> cogMap = functionMap.containsKey("COG") ? functionMap.get("COG") : null;
		HashMap<String, String> cogCateMap = functionMap.containsKey("COG category") ? functionMap.get("COG category")
				: null;

		if (cogMap != null && cogCateMap != null) {
			list.add("COG accession");
			list.add("COG name");
			list.add("COG category");
		}

		HashMap<String, String> nogMap = functionMap.containsKey("NOG") ? functionMap.get("NOG") : null;
		HashMap<String, String> nogCateMap = functionMap.containsKey("NOG category") ? functionMap.get("NOG category")
				: null;

		if (nogMap != null && nogCateMap != null) {
			list.add("NOG accession");
			list.add("NOG name");
			list.add("NOG category");
		}

		HashMap<String, String> keggMap = functionMap.containsKey("KEGG") ? functionMap.get("KEGG") : null;

		if (keggMap != null) {
			list.add("KEGG accession");
			list.add("KEGG name");
		}

		HashMap<String, String> koMap = functionMap.containsKey("KO") ? functionMap.get("KO") : null;
		if (koMap != null) {
			list.add("KO accession");
			list.add("KO name");
		}

		HashMap<String, String> goMap = functionMap.containsKey("GO") ? functionMap.get("GO") : null;
		if (goMap != null) {
			list.add("GOBP accession");
			list.add("GOBP name");

			list.add("GOCC accession");
			list.add("GOCC name");

			list.add("GOMF accession");
			list.add("GOMF name");
		}

		for (String title : list) {
			this.writer.println("<th>" + title + "</th>");
		}
		this.writer.println("</tr>");

		MetaProteinAnno1 delegrate = null;
		ArrayList<String> proNameList = new ArrayList<String>();

		for (MetaProteinAnno1 mp : proteins) {

			if (delegrate == null) {
				delegrate = mp;
				proNameList.add(mp.getPro().getName());
			} else {
				if (delegrate.getPro().getGroupId() == mp.getPro().getGroupId()) {
					proNameList.add(mp.getPro().getName());

				} else {

					this.writer.println("<tr>");

					this.writer.println("<td class=\"group-id\">" + (delegrate.getPro().getGroupId()) + "</td>");

					StringBuilder namesb = new StringBuilder("<td class=\"protein-names\">");
					for (String name : proNameList) {
						namesb.append(name).append(";").append("<br>");
					}

					namesb.append("</td>");
					this.writer.println(namesb);

					if (cogMap != null && cogCateMap != null) {
						String cog = delegrate.getCOG();
						String cogName = cogMap.get(cog);
						String cogCate = cogCateMap.get(cog);

						if (cog != null) {
							this.writer.println("<td><a href=\"" + getCOGLink(cog) + "\">" + cog + "</a></td>");
						} else {
							this.writer.println("<td class=\"empty-cell\">-</td>");
						}

						if (cogName != null) {
							this.writer.println("<td>" + cogName + "</td>");
						} else {
							this.writer.println("<td class=\"empty-cell\">-</td>");
						}

						if (cogCate != null) {
							this.writer.println("<td>" + cogCate + "</td>");
						} else {
							this.writer.println("<td class=\"empty-cell\">-</td>");
						}
					}

					if (nogMap != null && nogCateMap != null) {
						String nog = delegrate.getNOG();
						String nogName = nogMap.get(nog);
						String nogCate = nogCateMap.get(nog);

						if (nog != null) {
							this.writer.println("<td>" + nog + "</td>");
						} else {
							this.writer.println("<td class=\"empty-cell\">-</td>");
						}

						if (nogName != null) {
							this.writer.println("<td>" + nogName + "</td>");
						} else {
							this.writer.println("<td class=\"empty-cell\">-</td>");
						}

						if (nogCate != null) {
							this.writer.println("<td>" + nogCate + "</td>");
						} else {
							this.writer.println("<td class=\"empty-cell\">-</td>");
						}
					}

					if (keggMap != null) {
						String kegg = delegrate.getKEGG();
						String keggName = keggMap.get(kegg);

						if (kegg != null) {
							this.writer
									.println("<td><a href=\"" + getKEGGLink(kegg) + "\">" + kegg + "</a></td>");
						} else {
							this.writer.println("<td class=\"empty-cell\">-</td>");
						}

						if (keggName != null) {
							this.writer.println("<td>" + keggName + "</td>");
						} else {
							this.writer.println("<td class=\"empty-cell\">-</td>");
						}
					}

					if (koMap != null) {
						String ko = delegrate.getKO();
						String koName = koMap.get(ko);

						if (ko != null) {
							this.writer.println("<td><a href=\"" + getKEGGLink(ko) + "\">" + ko + "</a></td>");
						} else {
							this.writer.println("<td class=\"empty-cell\">-</td>");
						}

						if (koName != null) {
							this.writer.println("<td>" + koName + "</td>");
						} else {
							this.writer.println("<td class=\"empty-cell\">-</td>");
						}
					}

					if (goMap != null) {

						String[] gobp = delegrate.getGOBP();
						String[] gocc = delegrate.getGOCC();
						String[] gomf = delegrate.getGOMF();

						if (gobp != null && gobp.length > 0) {
							StringBuilder sb1 = new StringBuilder("<td><ul>");
							StringBuilder sb2 = new StringBuilder("<td><ul>");
							for (String go : gobp) {
								sb1.append("<li><a href=\"").append(getGOBPQuickGoLink(go)).append("\">").append(go)
										.append("</a></li>");
								if (goMap.containsKey(go)) {
									sb2.append("<li><a href=\"").append(getGOBPAmiGoLink(go)).append("\">")
											.append(goMap.get(go)).append("</a></li>");
								}
							}
							sb1.append("</ul></td>");
							sb2.append("</ul></td>");

							this.writer.println(sb1);
							this.writer.println(sb2);

						} else {
							this.writer.println("<td class=\"empty-cell\">-</td>");
							this.writer.println("<td class=\"empty-cell\">-</td>");
						}

						if (gocc != null && gocc.length > 0) {
							StringBuilder sb1 = new StringBuilder("<td><ul>");
							StringBuilder sb2 = new StringBuilder("<td><ul>");
							for (String go : gocc) {
								sb1.append("<li><a href=\"").append(getGOCCQuickGoLink(go)).append("\">").append(go)
										.append("</a></li>");
								if (goMap.containsKey(go)) {
									sb2.append("<li><a href=\"").append(getGOCCAmiGoLink(go)).append("\">")
											.append(goMap.get(go)).append("</a></li>");
								}
							}
							sb1.append("</ul></td>");
							sb2.append("</ul></td>");

							this.writer.println(sb1);
							this.writer.println(sb2);

						} else {
							this.writer.println("<td class=\"empty-cell\">-</td>");
							this.writer.println("<td class=\"empty-cell\">-</td>");
						}

						if (gomf != null && gomf.length > 0) {
							StringBuilder sb1 = new StringBuilder("<td><ul>");
							StringBuilder sb2 = new StringBuilder("<td><ul>");
							for (String go : gomf) {
								sb1.append("<li><a href=\"").append(getGOMFQuickGoLink(go)).append("\">").append(go)
										.append("</a></li>");
								if (goMap.containsKey(go)) {
									sb2.append("<li><a href=\"").append(getGOMFAmiGoLink(go)).append("\">")
											.append(goMap.get(go)).append("</a></li>");
								}
							}
							sb1.append("</ul></td>");
							sb2.append("</ul></td>");

							this.writer.println(sb1);
							this.writer.println(sb2);

						} else {
							this.writer.println("<td class=\"empty-cell\">-</td>");
							this.writer.println("<td class=\"empty-cell\">-</td>");
						}
					}

					this.writer.println("</tr>");

					delegrate = mp;
					proNameList = new ArrayList<String>();
					proNameList.add(mp.getPro().getName());
				}
			}
		}
	}

	public void close() {
		this.writer.write("\n</table>\n");
		this.writer.write("        </div>\n");  // Close table-container
		this.writer.write("    </div>\n");      // Close content-wrapper
		this.writer.write("    <div class=\"footer\">\n");
		this.writer.write("        <div class=\"logo\"><i class=\"fas fa-dna icon\"></i>MetaPilot</div>\n");
		this.writer.write("        <p>Advanced Metaproteomics Analysis Platform</p>\n");
		this.writer.write("        <p><i class=\"fas fa-university icon\"></i>Quadram Institute Bioscience</p>\n");
		LocalDateTime now = LocalDateTime.now();
		DateTimeFormatter formatter = DateTimeFormatter.ofPattern("yyyy-MM-dd HH:mm:ss");
		this.writer.write("        <p style=\"font-size: 0.9rem; margin-top: 10px;\">Generated on " + 
			now.format(formatter) + "</p>\n");
		this.writer.write("    </div>\n");      // Close footer
		this.writer.write("</div>\n");          // Close container
		this.writer.write("</body>\n</html>");
		this.writer.close();
	}
}
