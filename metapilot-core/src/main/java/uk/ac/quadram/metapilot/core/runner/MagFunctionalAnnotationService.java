package uk.ac.quadram.metapilot.core.runner;

import uk.ac.quadram.metapilot.core.function.v2.EnzymeCommission;
import uk.ac.quadram.metapilot.core.function.v2.GoObo;
import uk.ac.quadram.metapilot.task.io.pro.MetaProtein;
import uk.ac.quadram.metapilot.task.mag.MagDbItem;
import uk.ac.quadram.metapilot.task.mag.MagFuncSearcher;
import uk.ac.quadram.metapilot.task.mag.MetaProteinAnnoMag;

import java.io.File;
import java.io.IOException;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;

/**
 * Shared MAG functional annotation and static report generation.
 * <p>
 * The formal-release path keeps annotation data in memory and writes direct
 * TSV/HTML outputs. {@code functions.xml} is not part of this service contract.
 */
public final class MagFunctionalAnnotationService {

	private MagFunctionalAnnotationService() {
	}

	public static final class Request {
		public File resultDir;
		public MagDbItem magDb;
		public String microDbCatalogueId;
		public File functionDefinitionDb;
		public String analysisType;
		public String quantitativeMode;
		public String[] intensityTitles;
		public MetaProtein[] metaProteins;
		public String[] proteins;
		public HashMap<String, double[]> proteinIntensityMap;
		public HashMap<String, HashSet<String>> razorPepCountMap;
		public HashMap<String, Integer> razorPsmCountMap;
		/**
		 * Opt-in: also write the MetaAtlas-aligned function_abundance.parquet (see
		 * {@link FunctionalAnnotationExporter#exportMetaAtlasFunctionAbundance}). Defaults to
		 * false, so existing callers (in particular GUI tasks) get no behavior change unless
		 * they explicitly opt in. Intended to be set from a CLI flag (e.g.
		 * --export-metaatlas), not from GUI code.
		 */
		public boolean exportMetaAtlasParquet;
	}

	public static final class Result {
		public File functionalAnnotationDir;
		public File annotatedProteinsTsv;
		public File functionsHtml;
		public MetaProteinAnnoMag[] annotatedProteins;
	}

	public static Result run(Request request) throws SQLException, IOException {
		if (request == null) throw new IllegalArgumentException("request is null");
		if (request.resultDir == null) throw new IOException("resultDir is null");
		if (request.magDb == null) throw new IOException("magDb is null");
		if (request.functionDefinitionDb == null || !request.functionDefinitionDb.isFile()) {
			throw new IOException("Function definition database not found: "
					+ (request.functionDefinitionDb == null ? "null" : request.functionDefinitionDb.getAbsolutePath()));
		}
		File funcDir = new File(request.resultDir, "functional_annotation");
		funcDir.mkdirs();

		MetaProtein[] metaProteins = request.metaProteins != null ? request.metaProteins
				: buildMetaProteins(request.proteins, request.proteinIntensityMap,
						request.razorPepCountMap, request.razorPsmCountMap);
		String[] intensityTitles = request.intensityTitles != null ? request.intensityTitles : new String[0];
		FunctionalAnnotationResult annotation = annotateProteins(request.magDb, request.functionDefinitionDb,
				request.microDbCatalogueId, metaProteins);
		FunctionalAnnotationExporter.Result export = FunctionalAnnotationExporter.write(funcDir, annotation, intensityTitles);

		if (request.exportMetaAtlasParquet) {
			try {
				FunctionalAnnotationExporter.exportMetaAtlasFunctionAbundance(funcDir, annotation, intensityTitles);
			} catch (Exception e) {
				// Never let the new, guide-aligned export break the existing report generation.
			}
		}

		Result result = new Result();
		result.functionalAnnotationDir = funcDir;
		result.annotatedProteinsTsv = export.annotatedProteinsTsv;
		result.functionsHtml = export.functionsHtml;
		result.annotatedProteins = annotation.annotatedProteins;
		return result;
	}

	private static MetaProtein[] buildMetaProteins(String[] proteins, HashMap<String, double[]> proteinIntensityMap,
			HashMap<String, HashSet<String>> razorPepCountMap, HashMap<String, Integer> razorPsmCountMap) {
		if (proteins == null) return new MetaProtein[0];
		ArrayList<MetaProtein> list = new ArrayList<>();
		for (int i = 0; i < proteins.length; i++) {
			String protein = proteins[i];
			double[] intensity = proteinIntensityMap != null ? proteinIntensityMap.get(protein) : null;
			if (intensity == null) intensity = new double[0];
			int pepCount = razorPepCountMap != null && razorPepCountMap.containsKey(protein)
					? razorPepCountMap.get(protein).size() : 0;
			int psmCount = razorPsmCountMap != null && razorPsmCountMap.containsKey(protein)
					? razorPsmCountMap.get(protein) : 0;
			list.add(new MetaProtein(i, 0, protein, pepCount, psmCount, 0, intensity));
		}
		return list.toArray(new MetaProtein[0]);
	}

	private static FunctionalAnnotationResult annotateProteins(MagDbItem magDb, File funcDbFile, String microDbCatalogueId,
			MetaProtein[] proteins) throws SQLException, IOException {
		if (proteins == null) proteins = new MetaProtein[0];
		HashMap<String, ArrayList<MetaProtein>> genomeProMap = new HashMap<>();
		for (MetaProtein protein : proteins) {
			if (protein == null || protein.getName() == null || !protein.getName().startsWith(magDb.getIdentifier())) {
				continue;
			}
			String genome = protein.getName().split("_")[0];
			if (genome.endsWith(".1")) genome = genome.substring(0, genome.length() - 2);
			genomeProMap.computeIfAbsent(genome, k -> new ArrayList<>()).add(protein);
		}

		HashMap<String, String[]> usedCogMap = new HashMap<>();
		HashMap<String, String[]> usedNogMap = new HashMap<>();
		HashMap<String, String> usedKeggMap = new HashMap<>();
		HashMap<String, GoObo> usedGoMap = new HashMap<>();
		HashMap<String, EnzymeCommission> usedEcMap = new HashMap<>();
		ArrayList<MetaProteinAnnoMag> annotated = new ArrayList<>();

		for (String genome : genomeProMap.keySet()) {
			MetaProtein[] batch = genomeProMap.get(genome).toArray(new MetaProtein[0]);
			MagFuncSearcher searcher = new MagFuncSearcher(microDbCatalogueId, funcDbFile.getAbsolutePath(), magDb,
					genome);
			MetaProteinAnnoMag[] matches = searcher.match(batch);
			if (matches != null) {
				for (MetaProteinAnnoMag match : matches) annotated.add(match);
				usedCogMap.putAll(searcher.getUsedCogMap());
				usedNogMap.putAll(searcher.getUsedNogMap());
				usedKeggMap.putAll(searcher.getUsedKeggMap());
				usedEcMap.putAll(searcher.getUsedEcMap());
				usedGoMap.putAll(searcher.getUsedGoMap());
			}
		}

		FunctionalAnnotationResult result = new FunctionalAnnotationResult();
		result.annotatedProteins = annotated.toArray(new MetaProteinAnnoMag[0]);
		result.usedCogMap = usedCogMap;
		result.usedNogMap = usedNogMap;
		result.usedKeggMap = usedKeggMap;
		result.usedGoMap = usedGoMap;
		result.usedEcMap = usedEcMap;
		return result;
	}
}
