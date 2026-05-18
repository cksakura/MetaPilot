package uk.ac.quadram.metapilot.core.runner;

import uk.ac.quadram.metapilot.core.function.v2.EnzymeCommission;
import uk.ac.quadram.metapilot.core.function.v2.GoObo;
import uk.ac.quadram.metapilot.task.io.pro.MetaProtein;
import uk.ac.quadram.metapilot.task.io.pro.MetaProteinHtmlWriter;
import uk.ac.quadram.metapilot.task.mag.MetaProteinAnnoMag;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Arrays;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;

/**
 * Direct functional annotation report writer.
 */
public final class FunctionalAnnotationExporter {
	private FunctionalAnnotationExporter() {
	}

	public static final class Result {
		public File annotatedProteinsTsv;
		public File functionsHtml;
	}

	public static Result write(File funcDir, FunctionalAnnotationResult annotation, String[] intensityTitles)
			throws IOException {
		if (funcDir == null)
			throw new IOException("funcDir is null");
		if (annotation == null)
			annotation = new FunctionalAnnotationResult();
		if (annotation.annotatedProteins == null)
			annotation.annotatedProteins = new MetaProteinAnnoMag[0];
		if (intensityTitles == null)
			intensityTitles = new String[0];
		funcDir.mkdirs();

		File annotatedProteins = new File(funcDir, "Annotated_proteins.tsv");
		writeAnnotatedProteins(annotatedProteins, annotation, intensityTitles);

		writeGoReport(new File(funcDir, "Gene_Ontology_id_report.tsv"), annotation, intensityTitles);
		writeEcReport(new File(funcDir, "EC_id_report.tsv"), annotation, intensityTitles);
		writeKeggReport(new File(funcDir, "KEGG_ko_report.tsv"), annotation, intensityTitles);
		writeBriteReport(new File(funcDir, "BRITE_report.tsv"), annotation, intensityTitles);
		writeCogReport(new File(funcDir, "COG_accession_report.tsv"), annotation, intensityTitles);
		writeNogReport(new File(funcDir, "NOG_accession_report.tsv"), annotation, intensityTitles);
		writeCogCategoryReport(new File(funcDir, "COG_category_report.tsv"), annotation, intensityTitles);
		writeNogCategoryReport(new File(funcDir, "NOG_category_report.tsv"), annotation, intensityTitles);

		File html = new File(funcDir, "functions.html");
		MetaProteinHtmlWriter htmlWriter = new MetaProteinHtmlWriter(html);
		htmlWriter.write(annotation.annotatedProteins, annotation.usedCogMap, annotation.usedNogMap);
		htmlWriter.close();

		Result result = new Result();
		result.annotatedProteinsTsv = annotatedProteins;
		result.functionsHtml = html;
		return result;
	}

	private static void writeAnnotatedProteins(File out, FunctionalAnnotationResult annotation,
			String[] intensityTitles) throws IOException {
		MetaProteinAnnoMag[] proteins = annotation.annotatedProteins.clone();
		Arrays.sort(proteins, Comparator.comparingInt((MetaProteinAnnoMag p) -> p.getPro().getGroupId())
				.thenComparingInt(p -> p.getPro().getProteinId()).thenComparing(p -> p.getPro().getName()));
		try (PrintWriter pw = new PrintWriter(out)) {
			pw.print("Group_ID\tName\tPeptide count\tPSM count");
			for (String title : intensityTitles) {
				pw.print("\t");
				pw.print(normalizeIntensityTitle(title));
			}
			pw.println("\tProtein name\tSeed ortholog\tDescription\tTaxonomy Id\tTaxonomy name\tPreferred name"
					+ "\tGene_Ontology_id\tGene_Ontology_name\tGene_Ontology_namespace"
					+ "\tEC_id\tEC_description\tEC_reaction\tEC_synonyms"
					+ "\tKEGG_ko\tKEGG_Pathway_Entry\tKEGG_Pathway_Name\tKEGG_Module\tKEGG_Reaction\tKEGG_rclass"
					+ "\tBRITE\tKEGG_TC\tCAZy\tBiGG_Reaction\tPFAMs" + "\tCOG accession\tCOG category\tCOG name"
					+ "\tNOG accession\tNOG category\tNOG name");
			for (MetaProteinAnnoMag protein : proteins) {
				MetaProtein p = protein.getPro();
				StringBuilder sb = new StringBuilder();
				sb.append(p.getGroupId()).append("\t");
				sb.append(clean(p.getName())).append("\t");
				sb.append(p.getPepCount()).append("\t");
				sb.append(p.getPsmCount());
				for (double intensity : p.getIntensities()) {
					sb.append("\t").append(intensity);
				}
				for (String value : fixedFunctionalColumns(protein, annotation)) {
					sb.append("\t").append(clean(value));
				}
				pw.println(sb);
			}
		}
	}

	private static String[] fixedFunctionalColumns(MetaProteinAnnoMag protein, FunctionalAnnotationResult annotation) {
		return new String[] { protein.getPro_name(), protein.getSeed(), protein.getDescription(), protein.getTax_id(),
				protein.getBest_tax_level(), protein.getPrefered_name(), protein.getGene_Ontology(),
				joinGoNames(protein.getGene_Ontology(), annotation.usedGoMap),
				joinGoNamespaces(protein.getGene_Ontology(), annotation.usedGoMap), protein.getEC(),
				joinEcDescriptions(protein.getEC(), annotation.usedEcMap),
				joinEcReactions(protein.getEC(), annotation.usedEcMap),
				joinEcSynonyms(protein.getEC(), annotation.usedEcMap), protein.getKEGG_ko(), protein.getKEGG_Pathway(),
				joinMapValues(protein.getKEGG_Pathway(), annotation.usedKeggMap), protein.getKEGG_Module(),
				protein.getKEGG_Reaction(), protein.getKEGG_rclass(), protein.getBRITE(), protein.getKEGG_TC(),
				protein.getCAZy(), protein.getBiGG_Reaction(), protein.getPfam(), protein.getCog(),
				joinDefinitionColumn(protein.getCog(), annotation.usedCogMap, 0),
				joinDefinitionColumn(protein.getCog(), annotation.usedCogMap, 1), protein.getNog(),
				joinNogDefinitionColumn(protein.getNog(), annotation.usedNogMap, 0),
				joinNogDefinitionColumn(protein.getNog(), annotation.usedNogMap, 1) };
	}

	private static void writeGoReport(File out, FunctionalAnnotationResult annotation, String[] titles)
			throws IOException {
		HashMap<String, double[]> sums = aggregate(annotation.annotatedProteins, titles.length,
				p -> p.getGene_Ontology());
		try (PrintWriter pw = new PrintWriter(out)) {
			writeHeader(pw, "Gene_Ontology_id\tGene_Ontology_name\tGene_Ontology_namespace", titles);
			for (String id : sums.keySet()) {
				GoObo go = annotation.usedGoMap.get(id);
				pw.print(id + "\t" + clean(go != null ? go.getName() : "") + "\t"
						+ clean(go != null ? go.getNamespace() : ""));
				writeValues(pw, sums.get(id));
			}
		}
	}

	private static void writeEcReport(File out, FunctionalAnnotationResult annotation, String[] titles)
			throws IOException {
		HashMap<String, double[]> sums = aggregate(annotation.annotatedProteins, titles.length, p -> p.getEC());
		try (PrintWriter pw = new PrintWriter(out)) {
			writeHeader(pw, "EC_id\tEC_description\tEC_reaction\tEC_synonyms", titles);
			for (String id : sums.keySet()) {
				EnzymeCommission ec = annotation.usedEcMap.get(id);
				String synonyms = ec != null && ec.getAn() != null ? String.join(";", ec.getAn()) : "";
				pw.print(id + "\t" + clean(ec != null ? ec.getDe() : "") + "\t" + clean(ec != null ? ec.getCa() : "")
						+ "\t" + clean(synonyms));
				writeValues(pw, sums.get(id));
			}
		}
	}

	private static void writeKeggReport(File out, FunctionalAnnotationResult annotation, String[] titles)
			throws IOException {
		HashMap<String, double[]> sums = aggregate(annotation.annotatedProteins, titles.length, p -> p.getKEGG_ko());
		try (PrintWriter pw = new PrintWriter(out)) {
			writeHeader(pw, "KEGG_ko\tKEGG_ko_description", titles);
			for (String id : sums.keySet()) {
				pw.print(id + "\t" + clean(annotation.usedKeggMap.get(id)));
				writeValues(pw, sums.get(id));
			}
		}
	}

	private static void writeBriteReport(File out, FunctionalAnnotationResult annotation, String[] titles)
			throws IOException {
		writeSimpleReport(out, "BRITE", aggregate(annotation.annotatedProteins, titles.length, p -> p.getBRITE()),
				titles);
	}

	private static void writeCogReport(File out, FunctionalAnnotationResult annotation, String[] titles)
			throws IOException {
		HashMap<String, double[]> sums = aggregate(annotation.annotatedProteins, titles.length, p -> p.getCog());
		try (PrintWriter pw = new PrintWriter(out)) {
			writeHeader(pw, "COG accession\tCOG category\tCOG name", titles);
			for (String id : sums.keySet()) {
				String[] def = annotation.usedCogMap.get(id);
				pw.print(id + "\t" + clean(def != null && def.length > 0 ? def[0] : "") + "\t"
						+ clean(def != null && def.length > 1 ? def[1] : ""));
				writeValues(pw, sums.get(id));
			}
		}
	}

	private static void writeNogReport(File out, FunctionalAnnotationResult annotation, String[] titles)
			throws IOException {
		HashMap<String, double[]> sums = aggregate(annotation.annotatedProteins, titles.length, p -> p.getNog());
		try (PrintWriter pw = new PrintWriter(out)) {
			writeHeader(pw, "NOG accession\tNOG category\tNOG name", titles);
			for (String id : sums.keySet()) {
				String lookupId = id.contains("|") ? id.substring(0, id.indexOf("|")) : id;
				String[] def = annotation.usedNogMap.get(lookupId);
				pw.print(id + "\t" + clean(def != null && def.length > 0 ? def[0] : "") + "\t"
						+ clean(def != null && def.length > 1 ? def[1] : ""));
				writeValues(pw, sums.get(id));
			}
		}
	}

	private static void writeCogCategoryReport(File out, FunctionalAnnotationResult annotation, String[] titles)
			throws IOException {
		HashMap<String, double[]> sums = new HashMap<>();
		for (MetaProteinAnnoMag protein : annotation.annotatedProteins) {
			String cog = clean(protein.getCog());
			String[] def = annotation.usedCogMap.get(cog);
			addSplitIntensities(sums, def != null && def.length > 0 ? def[0] : "", protein.getPro().getIntensities(),
					titles.length);
		}
		writeSimpleReport(out, "COG Category", sums, titles);
	}

	private static void writeNogCategoryReport(File out, FunctionalAnnotationResult annotation, String[] titles)
			throws IOException {
		HashMap<String, double[]> sums = new HashMap<>();
		for (MetaProteinAnnoMag protein : annotation.annotatedProteins) {
			for (String nog : splitIds(protein.getNog())) {
				String lookupId = nog.contains("|") ? nog.substring(0, nog.indexOf("|")) : nog;
				String[] def = annotation.usedNogMap.get(lookupId);
				addSplitIntensities(sums, def != null && def.length > 0 ? def[0] : "",
						protein.getPro().getIntensities(), titles.length);
			}
		}
		writeSimpleReport(out, "NOG Category", sums, titles);
	}

	private static void writeSimpleReport(File out, String firstColumn, HashMap<String, double[]> sums, String[] titles)
			throws IOException {
		try (PrintWriter pw = new PrintWriter(out)) {
			writeHeader(pw, firstColumn, titles);
			for (String id : sums.keySet()) {
				pw.print(id);
				writeValues(pw, sums.get(id));
			}
		}
	}

	private static HashMap<String, double[]> aggregate(MetaProteinAnnoMag[] proteins, int width,
			IdExtractor extractor) {
		HashMap<String, double[]> sums = new HashMap<>();
		for (MetaProteinAnnoMag protein : proteins) {
			String[] ids = splitIds(extractor.ids(protein));
			if (ids.length == 0)
				continue;
			double[] intensity = protein.getPro().getIntensities();
			double share = 1.0 / ids.length;
			for (String id : ids) {
				double[] target = sums.computeIfAbsent(id, k -> new double[width]);
				for (int i = 0; i < Math.min(width, intensity.length); i++) {
					target[i] += intensity[i] * share;
				}
			}
		}
		return sums;
	}

	private static void addSplitIntensities(Map<String, double[]> sums, String ids, double[] intensity, int width) {
		String[] split = splitIds(ids);
		if (split.length == 0)
			return;
		double share = 1.0 / split.length;
		for (String id : split) {
			double[] target = sums.computeIfAbsent(id, k -> new double[width]);
			for (int i = 0; i < Math.min(width, intensity.length); i++) {
				target[i] += intensity[i] * share;
			}
		}
	}

	private static String[] splitIds(String ids) {
		if (ids == null || ids.trim().isEmpty() || "-".equals(ids.trim()))
			return new String[0];
		return Arrays.stream(ids.split(",")).map(String::trim).filter(s -> !s.isEmpty() && !"-".equals(s))
				.toArray(String[]::new);
	}

	private static String joinGoNames(String ids, HashMap<String, GoObo> map) {
		StringBuilder sb = new StringBuilder();
		for (String id : splitIds(ids)) {
			GoObo go = map.get(id);
			appendDelimited(sb, go != null ? go.getName() : "");
		}
		return sb.toString();
	}

	private static String joinGoNamespaces(String ids, HashMap<String, GoObo> map) {
		StringBuilder sb = new StringBuilder();
		for (String id : splitIds(ids)) {
			GoObo go = map.get(id);
			appendDelimited(sb, go != null ? go.getNamespace() : "");
		}
		return sb.toString();
	}

	private static String joinEcDescriptions(String ids, HashMap<String, EnzymeCommission> map) {
		StringBuilder sb = new StringBuilder();
		for (String id : splitIds(ids)) {
			EnzymeCommission ec = map.get(id);
			appendDelimited(sb, ec != null ? ec.getDe() : "");
		}
		return sb.toString();
	}

	private static String joinEcReactions(String ids, HashMap<String, EnzymeCommission> map) {
		StringBuilder sb = new StringBuilder();
		for (String id : splitIds(ids)) {
			EnzymeCommission ec = map.get(id);
			appendDelimited(sb, ec != null ? ec.getCa() : "");
		}
		return sb.toString();
	}

	private static String joinEcSynonyms(String ids, HashMap<String, EnzymeCommission> map) {
		StringBuilder sb = new StringBuilder();
		HashSet<String> seen = new HashSet<>();
		for (String id : splitIds(ids)) {
			EnzymeCommission ec = map.get(id);
			if (ec == null || ec.getAn() == null)
				continue;
			for (String synonym : ec.getAn()) {
				if (seen.add(synonym))
					appendDelimited(sb, synonym);
			}
		}
		return sb.toString();
	}

	private static String joinMapValues(String ids, HashMap<String, String> map) {
		StringBuilder sb = new StringBuilder();
		for (String id : splitIds(ids)) {
			appendDelimited(sb, map.get(id));
		}
		return sb.toString();
	}

	private static String joinDefinitionColumn(String ids, HashMap<String, String[]> map, int index) {
		StringBuilder sb = new StringBuilder();
		for (String id : splitIds(ids)) {
			String[] def = map.get(id);
			appendDelimited(sb, def != null && def.length > index ? def[index] : "");
		}
		return sb.toString();
	}

	private static String joinNogDefinitionColumn(String ids, HashMap<String, String[]> map, int index) {
		StringBuilder sb = new StringBuilder();
		for (String id : splitIds(ids)) {
			String lookupId = id.contains("|") ? id.substring(0, id.indexOf("|")) : id;
			String[] def = map.get(lookupId);
			appendDelimited(sb, def != null && def.length > index ? def[index] : "");
		}
		return sb.toString();
	}

	private static void appendDelimited(StringBuilder sb, String value) {
		if (value == null || value.isEmpty())
			return;
		if (sb.length() > 0)
			sb.append(";");
		sb.append(value);
	}

	private static void writeHeader(PrintWriter pw, String firstColumns, String[] titles) {
		pw.print(firstColumns);
		for (String title : titles) {
			pw.print("\t");
			pw.print(normalizeIntensityTitle(title));
		}
		pw.println();
	}

	private static void writeValues(PrintWriter pw, double[] values) {
		for (double value : values) {
			pw.print("\t");
			pw.print(value);
		}
		pw.println();
	}

	private static String normalizeIntensityTitle(String title) {
		if (title == null)
			return "";
		return title.startsWith("Intensity") ? title : "Intensity " + title;
	}

	private static String clean(String value) {
		return value == null ? "" : value.replace('\t', ' ').replace('\n', ' ').replace('\r', ' ');
	}

	@FunctionalInterface
	private interface IdExtractor {
		String ids(MetaProteinAnnoMag protein);
	}
}
