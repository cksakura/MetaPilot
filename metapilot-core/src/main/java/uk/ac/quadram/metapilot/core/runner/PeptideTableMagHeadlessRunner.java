package uk.ac.quadram.metapilot.core.runner;

import uk.ac.quadram.metapilot.core.dia.MagPredictedPeptideProteinMapper;
import uk.ac.quadram.metapilot.core.mag.DdaMagPostAnalysisService;
import uk.ac.quadram.metapilot.core.workflow.ProgressCallback;
import uk.ac.quadram.metapilot.task.io.pep.MetaPeptide;
import uk.ac.quadram.metapilot.task.mag.MagDbConfigIO;
import uk.ac.quadram.metapilot.task.mag.MagDbItem;
import uk.ac.quadram.metapilot.task.mag.MagDbItemResolver;
import uk.ac.quadram.metapilot.task.par.MetaConstants;

import org.dom4j.Element;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.nio.charset.StandardCharsets;
import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.Locale;
import java.util.Map;
import java.util.Set;
import java.util.TreeSet;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * Headless MAG post-analysis starting from a plain peptide intensity table.
 *
 * <p>Input: a TSV/CSV table with one peptide column (header {@code peptide} or
 * {@code sequence}, case-insensitive) and intensity columns. Every other column whose
 * non-missing values are all numeric is taken as a sample intensity column; any other
 * column is ignored and reported in the log.
 *
 * <p>Steps:
 * <ol>
 * <li>Normalise sequences to plain residues (modification annotations, flanking
 * residues and terminal markers removed) and sum intensities of rows that collapse
 * to the same sequence.</li>
 * <li>Map sequences to MAG proteins by exact match against
 * {@code <MAG db>/predicted/*.predicted.tsv}
 * ({@link MagPredictedPeptideProteinMapper#mapPeptideProteins}).</li>
 * <li>Genome filter: a genome is accepted when it has at least
 * {@code minUniquePeptides} genome-unique peptides (peptides whose every mapped
 * protein belongs to that genome). The default 0 accepts every genome with at least
 * one mapped peptide. There is no score and no decoy in this input, so no
 * FDR-controlled genome inference is possible here.</li>
 * <li>Razor protein assignment and genome/protein intensity roll-up
 * ({@link DdaMagPostAnalysisService#analyzeWithAcceptedGenomes}), taxonomy
 * ({@link MagTaxonomyAnalysisService}) and functional annotation
 * ({@link MagFunctionalAnnotationService}) - the same services the DDA MAG runners
 * use.</li>
 * </ol>
 *
 * <p>Outputs under {@code outDir}: {@code final_peptides.tsv},
 * {@code razor_proteins.tsv}, {@code all_proteins.tsv}, {@code genome_filter.tsv},
 * {@code run_summary.tsv}, {@code taxonomy_analysis/} and (when a function database
 * is given) {@code functional_annotation/}.
 */
public final class PeptideTableMagHeadlessRunner {

	public static final String TASK_NAME = "peptide-table-mag";

	private static final Set<String> SEQUENCE_HEADERS = Set.of("peptide", "sequence");
	private static final Set<String> MISSING_VALUES = Set.of("", "na", "nan", "null", "none", "-", "#n/a");
	private static final Pattern FLANKED = Pattern.compile("^[A-Za-z\\-]\\.(.+)\\.[A-Za-z\\-]$");
	private static final Pattern BRACKETED = Pattern.compile("\\([^)]*\\)|\\[[^\\]]*\\]|\\{[^}]*\\}");

	private PeptideTableMagHeadlessRunner() {
	}

	public static final class Request {
		public File peptideTable;
		public File outDir;
		/** MAG registry JSON; null resolves through {@link MagDbItemResolver#resolveConfigFile}. */
		public File magdbConfigJson;
		public String catalogueId;
		/** Optional catalogue version; null uses the registry's used/current version. */
		public String catalogueVersion;
		/** Optional func_def.db; functional annotation is skipped when null or missing. */
		public File funcDbFile;
		/** Minimum genome-unique peptides for a genome to be accepted. 0 = no filter. */
		public int minUniquePeptides = 0;
		public int threads = 4;
		public boolean exportMetaAtlasParquet;
		public ProgressCallback progressCallback;
	}

	public static final class Result {
		public File outDir;
		public int inputRows;
		public int distinctSequences;
		public int mappedSequences;
		public int genomesWithMappedPeptides;
		public int acceptedGenomes;
		public int peptidesWithRazorProtein;
		public String[] sampleTitles;
		public DdaMagPostAnalysisService.Result analysis;
	}

	// ---------------------------------------------------------------- run

	public static Result run(Request request) throws Exception {
		if (request == null)
			throw new IllegalArgumentException("request is null");
		if (request.peptideTable == null || !request.peptideTable.isFile())
			throw new IOException("Peptide table not found: " + request.peptideTable);
		if (request.outDir == null)
			throw new IllegalArgumentException("outDir is null");
		if (request.minUniquePeptides < 0)
			throw new IllegalArgumentException("minUniquePeptides must be >= 0, got " + request.minUniquePeptides);
		ProgressCallback cb = request.progressCallback != null ? request.progressCallback
				: new ProgressCallback.NoOp();

		File cfg = MagDbItemResolver.resolveConfigFile(request.magdbConfigJson);
		MagDbItem mag = resolveMagDb(request.catalogueId, request.catalogueVersion, cfg);
		if (mag == null) {
			throw new IllegalArgumentException("MAG database not found in " + cfg + " for catalogue="
					+ request.catalogueId + " version=" + request.catalogueVersion);
		}
		File predicted = mag.getPredictedPepFile();
		if (predicted == null || !predicted.isDirectory()) {
			throw new IOException("MAG predicted peptide folder not found: " + predicted);
		}
		String magIdentifier = mag.getIdentifier();
		cb.onMessage(TASK_NAME + ": MAG DB " + mag.getCatalogueID() + " / "
				+ (mag.getUsedVersion() != null ? mag.getUsedVersion().getCatalogueVersion() : "?") + ", predicted="
				+ predicted.getAbsolutePath());

		File outDir = request.outDir;
		outDir.mkdirs();

		// 1. read table
		PeptideTable table = PeptideTable.read(request.peptideTable, cb);
		int fileCount = table.sampleTitles.length;
		cb.onMessage(TASK_NAME + ": " + table.inputRows + " rows -> " + table.intensityBySequence.size()
				+ " distinct sequences, " + fileCount + " intensity columns: " + String.join(", ", table.sampleTitles));

		// 2. map to MAG proteins
		cb.onMessage(TASK_NAME + ": mapping peptides to MAG proteins ...");
		Map<String, Set<String>> pepProteins = MagPredictedPeptideProteinMapper.mapPeptideProteins(
				table.intensityBySequence.keySet(), predicted, null, Math.max(1, request.threads));

		// 3. genome-unique counting and filter
		HashMap<String, HashSet<String>> genomeMapped = new HashMap<>();
		HashMap<String, HashSet<String>> genomeUnique = new HashMap<>();
		HashMap<String, String> uniqueGenomeOf = new HashMap<>();
		int mapped = 0;
		for (String seq : table.intensityBySequence.keySet()) {
			Set<String> proteins = pepProteins.get(seq);
			if (proteins == null || proteins.isEmpty())
				continue;
			mapped++;
			Set<String> genomes = genomesOf(proteins, magIdentifier);
			for (String g : genomes)
				genomeMapped.computeIfAbsent(g, k -> new HashSet<>()).add(seq);
			if (genomes.size() == 1 && allMag(proteins, magIdentifier)) {
				String g = genomes.iterator().next();
				genomeUnique.computeIfAbsent(g, k -> new HashSet<>()).add(seq);
				uniqueGenomeOf.put(seq, g);
			}
		}
		TreeSet<String> accepted = new TreeSet<>();
		for (String g : genomeMapped.keySet()) {
			int unique = genomeUnique.containsKey(g) ? genomeUnique.get(g).size() : 0;
			if (unique >= request.minUniquePeptides)
				accepted.add(g);
		}
		cb.onMessage(TASK_NAME + ": mapped sequences=" + mapped + "/" + table.intensityBySequence.size()
				+ ", genomes with mapped peptides=" + genomeMapped.size() + ", accepted (>= "
				+ request.minUniquePeptides + " genome-unique peptides)=" + accepted.size());
		writeGenomeFilter(new File(outDir, "genome_filter.tsv"), genomeMapped, genomeUnique, accepted);
		if (accepted.isEmpty()) {
			throw new IOException(TASK_NAME + ": no genome passed the filter (min genome-unique peptides = "
					+ request.minUniquePeptides + "); see genome_filter.tsv");
		}

		// 4. razor assignment + roll-up
		ArrayList<TablePeptide> peptideList = new ArrayList<>();
		for (Map.Entry<String, double[]> e : table.intensityBySequence.entrySet()) {
			Set<String> proteins = pepProteins.get(e.getKey());
			String[] proArray = proteins == null ? new String[0] : proteins.toArray(new String[0]);
			Arrays.sort(proArray);
			peptideList.add(new TablePeptide(e.getKey(), proArray, e.getValue()));
		}
		TablePeptide[] peptides = peptideList.toArray(new TablePeptide[0]);
		DdaMagPostAnalysisService.Result analysis = DdaMagPostAnalysisService.analyzeWithAcceptedGenomes(peptides,
				fileCount, magIdentifier, accepted);

		int withRazor = 0;
		for (int i = 0; i < peptides.length; i++) {
			if (i < analysis.razorProteins.length && analysis.razorProteins[i] != null)
				withRazor++;
		}

		// 5. tables
		writePeptides(new File(outDir, "final_peptides.tsv"), peptides, analysis.razorProteins, uniqueGenomeOf,
				table.inputRowCount, magIdentifier, table.sampleTitles);
		writeProteins(new File(outDir, "razor_proteins.tsv"), analysis.proteins, analysis.proteinIntensityMap,
				analysis.proteinTotalIntenMap, analysis.totalPepCountMap, analysis.razorPepCountMap, magIdentifier,
				table.sampleTitles);
		writeProteins(new File(outDir, "all_proteins.tsv"), analysis.allProteins, analysis.allProteinIntensityMap,
				analysis.allProteinTotalIntenMap, analysis.allProteinPepCountMap, analysis.razorPepCountMap,
				magIdentifier, table.sampleTitles);

		// 6. taxonomy
		if (analysis.genomeIntensityMap.isEmpty()) {
			cb.onMessage(TASK_NAME + ": taxonomy skipped - no genome received intensity");
		} else {
			cb.onMessage(TASK_NAME + ": taxonomy report export ...");
			MagTaxonomyAnalysisService.Request tax = new MagTaxonomyAnalysisService.Request();
			tax.taxonomyAnalysisDir = new File(outDir, "taxonomy_analysis");
			tax.magDb = mag;
			tax.generateFinalGenomeRanking = true;
			tax.genomeIntensityMap = analysis.genomeIntensityMap;
			tax.genomeScoreMap = null; // no genome score exists for this input
			tax.genomeStats = null;
			tax.psmCountMap = analysis.psmCountMap;
			tax.totalPepCountMap = analysis.totalPepCountMap;
			tax.razorPepCountMap = analysis.razorPepCountMap;
			tax.uniquePepCountMap = analysis.uniquePepCountMap;
			tax.totalPepWithILMap = analysis.totalPepWithILMap;
			tax.razorPepWithILMap = analysis.razorPepWithILMap;
			tax.uniquePepWithILMap = analysis.uniquePepWithILMap;
			tax.intensityTitles = table.sampleTitles;
			tax.progressCallback = cb;
			tax.exportMetaAtlasParquet = request.exportMetaAtlasParquet;
			MagTaxonomyAnalysisService.run(tax);
		}

		// 7. functions
		if (request.funcDbFile != null && request.funcDbFile.isFile()) {
			cb.onMessage(TASK_NAME + ": functional annotation export ...");
			MagFunctionalAnnotationService.Request fn = new MagFunctionalAnnotationService.Request();
			fn.resultDir = outDir;
			fn.magDb = mag;
			fn.microDbCatalogueId = mag.getCatalogueID();
			fn.functionDefinitionDb = request.funcDbFile;
			fn.analysisType = MetaConstants.pFindMag;
			fn.quantitativeMode = MetaConstants.labelFree;
			fn.intensityTitles = table.sampleTitles;
			fn.metaProteins = analysis.metaProteins;
			fn.exportMetaAtlasParquet = request.exportMetaAtlasParquet;
			MagFunctionalAnnotationService.run(fn);
		} else {
			cb.onMessage(TASK_NAME + ": functional annotation skipped (no function database given)");
		}

		Result result = new Result();
		result.outDir = outDir;
		result.inputRows = table.inputRows;
		result.distinctSequences = table.intensityBySequence.size();
		result.mappedSequences = mapped;
		result.genomesWithMappedPeptides = genomeMapped.size();
		result.acceptedGenomes = accepted.size();
		result.peptidesWithRazorProtein = withRazor;
		result.sampleTitles = table.sampleTitles;
		result.analysis = analysis;
		writeSummary(new File(outDir, "run_summary.tsv"), request, mag, table, result);

		cb.onMessage(TASK_NAME + ": peptides with razor protein=" + withRazor + "/" + peptides.length
				+ ", razor proteins=" + analysis.proteins.length + ", genomes with intensity="
				+ analysis.genomeIntensityMap.size());
		cb.onComplete(outDir.getAbsolutePath());
		return result;
	}

	// ---------------------------------------------------------------- input

	static final class PeptideTable {
		String[] sampleTitles;
		int inputRows;
		/** normalised sequence -> summed intensity per sample; insertion order = first appearance */
		LinkedHashMap<String, double[]> intensityBySequence = new LinkedHashMap<>();
		HashMap<String, Integer> inputRowCount = new HashMap<>();

		static PeptideTable read(File file, ProgressCallback cb) throws IOException {
			char delim = file.getName().toLowerCase(Locale.ROOT).endsWith(".csv") ? ',' : '\t';
			ArrayList<String[]> rows = new ArrayList<>();
			String[] header;
			try (BufferedReader reader = new BufferedReader(
					new InputStreamReader(new FileInputStream(file), StandardCharsets.UTF_8))) {
				String line = reader.readLine();
				if (line == null)
					throw new IOException("Peptide table is empty: " + file);
				if (!line.isEmpty() && line.charAt(0) == '﻿')
					line = line.substring(1);
				header = split(line, delim);
				while ((line = reader.readLine()) != null) {
					if (!line.isBlank())
						rows.add(split(line, delim));
				}
			}

			int seqCol = -1;
			for (int i = 0; i < header.length; i++) {
				if (SEQUENCE_HEADERS.contains(header[i].trim().toLowerCase(Locale.ROOT))) {
					if (seqCol >= 0)
						throw new IOException("Both '" + header[seqCol] + "' and '" + header[i]
								+ "' look like the peptide column; keep only one");
					seqCol = i;
				}
			}
			if (seqCol < 0)
				throw new IOException("No 'peptide' or 'sequence' column in " + file.getName() + " (header: "
						+ String.join(", ", header) + ")");

			ArrayList<Integer> intenCols = new ArrayList<>();
			ArrayList<String> ignored = new ArrayList<>();
			for (int c = 0; c < header.length; c++) {
				if (c == seqCol)
					continue;
				boolean numeric = true;
				for (String[] row : rows) {
					String v = c < row.length ? row[c].trim() : "";
					if (!isMissing(v) && parse(v) == null) {
						numeric = false;
						break;
					}
				}
				if (numeric)
					intenCols.add(c);
				else
					ignored.add(header[c]);
			}
			if (intenCols.isEmpty())
				throw new IOException("No numeric intensity column found in " + file.getName());
			if (!ignored.isEmpty())
				cb.onMessage(TASK_NAME + ": ignored non-numeric column(s): " + String.join(", ", ignored));

			PeptideTable t = new PeptideTable();
			t.sampleTitles = new String[intenCols.size()];
			for (int i = 0; i < intenCols.size(); i++)
				t.sampleTitles[i] = sampleTitle(header[intenCols.get(i)]);
			if (new HashSet<>(Arrays.asList(t.sampleTitles)).size() != t.sampleTitles.length)
				throw new IOException("Duplicate sample names after removing the 'Intensity' prefix: "
						+ String.join(", ", t.sampleTitles));

			int skipped = 0;
			for (String[] row : rows) {
				String raw = seqCol < row.length ? row[seqCol] : "";
				String seq = normaliseSequence(raw);
				if (seq.isEmpty()) {
					skipped++;
					continue;
				}
				t.inputRows++;
				double[] target = t.intensityBySequence.computeIfAbsent(seq, k -> new double[intenCols.size()]);
				for (int i = 0; i < intenCols.size(); i++) {
					int c = intenCols.get(i);
					Double v = c < row.length ? parse(row[c].trim()) : null;
					if (v != null && v > 0 && !v.isInfinite())
						target[i] += v;
				}
				t.inputRowCount.merge(seq, 1, Integer::sum);
			}
			if (skipped > 0)
				cb.onMessage(TASK_NAME + ": skipped " + skipped + " row(s) with an empty peptide sequence");
			return t;
		}
	}

	/**
	 * Plain residue sequence: flanking residues ({@code K.PEPTIDE.R}) and anything in
	 * (), [] or {} removed, then only upper-case letters kept. Lower-case letters are
	 * dropped because they are terminal markers in common formats ({@code n[43]PEPTIDE}).
	 */
	static String normaliseSequence(String raw) {
		if (raw == null)
			return "";
		String s = raw.trim();
		Matcher m = FLANKED.matcher(s);
		if (m.matches())
			s = m.group(1);
		s = BRACKETED.matcher(s).replaceAll("");
		StringBuilder sb = new StringBuilder(s.length());
		for (int i = 0; i < s.length(); i++) {
			char ch = s.charAt(i);
			if (ch >= 'A' && ch <= 'Z')
				sb.append(ch);
		}
		return sb.toString();
	}

	static String sampleTitle(String header) {
		String h = header.trim();
		String lower = h.toLowerCase(Locale.ROOT);
		if (lower.startsWith("intensity") && h.length() > "intensity".length()) {
			char sep = h.charAt("intensity".length());
			if (sep == ' ' || sep == '_' || sep == '.' || sep == '-' || sep == ':') {
				String rest = h.substring("intensity".length() + 1).trim();
				if (!rest.isEmpty())
					return rest;
			}
		}
		return h;
	}

	private static boolean isMissing(String v) {
		return MISSING_VALUES.contains(v.toLowerCase(Locale.ROOT));
	}

	private static Double parse(String v) {
		if (isMissing(v))
			return null;
		try {
			return Double.parseDouble(v);
		} catch (NumberFormatException e) {
			return null;
		}
	}

	/** Splits one line; double-quoted fields may contain the delimiter ("" is an escaped quote). */
	static String[] split(String line, char delim) {
		ArrayList<String> out = new ArrayList<>();
		StringBuilder sb = new StringBuilder();
		boolean quoted = false;
		for (int i = 0; i < line.length(); i++) {
			char ch = line.charAt(i);
			if (quoted) {
				if (ch == '"') {
					if (i + 1 < line.length() && line.charAt(i + 1) == '"') {
						sb.append('"');
						i++;
					} else {
						quoted = false;
					}
				} else {
					sb.append(ch);
				}
			} else if (ch == '"' && sb.length() == 0) {
				quoted = true;
			} else if (ch == delim) {
				out.add(sb.toString());
				sb.setLength(0);
			} else {
				sb.append(ch);
			}
		}
		out.add(sb.toString());
		return out.toArray(new String[0]);
	}

	// ---------------------------------------------------------------- output

	private static void writeGenomeFilter(File out, Map<String, HashSet<String>> mapped,
			Map<String, HashSet<String>> unique, Set<String> accepted) throws IOException {
		String[] genomes = mapped.keySet().toArray(new String[0]);
		Arrays.sort(genomes, (a, b) -> {
			int ua = unique.containsKey(a) ? unique.get(a).size() : 0;
			int ub = unique.containsKey(b) ? unique.get(b).size() : 0;
			if (ua != ub)
				return Integer.compare(ub, ua);
			int c = Integer.compare(mapped.get(b).size(), mapped.get(a).size());
			return c != 0 ? c : a.compareTo(b);
		});
		try (PrintWriter pw = new PrintWriter(out, StandardCharsets.UTF_8)) {
			pw.println("Genome\tMapped peptides\tGenome-unique peptides\tAccepted");
			for (String g : genomes) {
				int u = unique.containsKey(g) ? unique.get(g).size() : 0;
				pw.println(g + "\t" + mapped.get(g).size() + "\t" + u + "\t" + accepted.contains(g));
			}
		}
	}

	private static void writePeptides(File out, TablePeptide[] peptides, String[] razorProteins,
			Map<String, String> uniqueGenomeOf, Map<String, Integer> inputRowCount, String magIdentifier,
			String[] samples) throws IOException {
		DecimalFormat df = intensityFormat();
		try (PrintWriter pw = new PrintWriter(out, StandardCharsets.UTF_8)) {
			StringBuilder h = new StringBuilder(
					"Sequence\tInput rows\tRazor protein\tRazor genome\tUnique genome\tMapped genome count\t"
							+ "Mapped protein count\tMapped genomes\tMapped proteins\tIntensity");
			for (String s : samples)
				h.append("\tIntensity ").append(s);
			pw.println(h);
			for (int i = 0; i < peptides.length; i++) {
				TablePeptide p = peptides[i];
				String razor = i < razorProteins.length && razorProteins[i] != null ? razorProteins[i] : "";
				Set<String> genomes = genomesOf(Arrays.asList(p.getProteins()), magIdentifier);
				StringBuilder sb = new StringBuilder();
				sb.append(p.getSequence()).append('\t');
				sb.append(inputRowCount.getOrDefault(p.getSequence(), 1)).append('\t');
				sb.append(razor).append('\t');
				sb.append(razor.isEmpty() ? "" : genomeOf(razor)).append('\t');
				sb.append(uniqueGenomeOf.getOrDefault(p.getSequence(), "")).append('\t');
				sb.append(genomes.size()).append('\t');
				sb.append(p.getProteins().length).append('\t');
				sb.append(String.join(";", genomes)).append('\t');
				sb.append(String.join(";", p.getProteins())).append('\t');
				double[] inten = p.getIntensity();
				sb.append(df.format(sum(inten)));
				for (double v : inten)
					sb.append('\t').append(df.format(v));
				pw.println(sb);
			}
		}
	}

	private static void writeProteins(File out, String[] proteins, Map<String, double[]> intensityMap,
			Map<String, Double> totalMap, Map<String, HashSet<String>> allPepMap,
			Map<String, HashSet<String>> razorPepMap, String magIdentifier, String[] samples) throws IOException {
		DecimalFormat df = intensityFormat();
		try (PrintWriter pw = new PrintWriter(out, StandardCharsets.UTF_8)) {
			StringBuilder h = new StringBuilder(
					"Protein\tGenome\tPeptide count (all)\tPeptide count (razor)\tIntensity");
			for (String s : samples)
				h.append("\tIntensity ").append(s);
			h.append("\tPeptides (all)\tPeptides (razor)");
			pw.println(h);
			for (String protein : proteins) {
				HashSet<String> all = allPepMap.getOrDefault(protein, new HashSet<>());
				HashSet<String> razor = razorPepMap.getOrDefault(protein, new HashSet<>());
				StringBuilder sb = new StringBuilder();
				sb.append(protein).append('\t');
				sb.append(isMag(protein, magIdentifier) ? genomeOf(protein) : "").append('\t');
				sb.append(all.size()).append('\t');
				sb.append(razor.size()).append('\t');
				sb.append(df.format(totalMap.getOrDefault(protein, 0.0)));
				double[] inten = intensityMap.getOrDefault(protein, new double[0]);
				for (int j = 0; j < samples.length; j++)
					sb.append('\t').append(df.format(j < inten.length ? inten[j] : 0.0));
				sb.append('\t').append(String.join(";", new TreeSet<>(all)));
				sb.append('\t').append(String.join(";", new TreeSet<>(razor)));
				pw.println(sb);
			}
		}
	}

	private static void writeSummary(File out, Request request, MagDbItem mag, PeptideTable table, Result r)
			throws IOException {
		try (PrintWriter pw = new PrintWriter(out, StandardCharsets.UTF_8)) {
			pw.println("Item\tValue");
			pw.println("Input table\t" + request.peptideTable.getAbsolutePath());
			pw.println("MAG catalogue\t" + mag.getCatalogueID());
			pw.println("MAG version\t" + (mag.getUsedVersion() != null ? mag.getUsedVersion().getCatalogueVersion() : ""));
			pw.println("Samples\t" + String.join(";", table.sampleTitles));
			pw.println("Min genome-unique peptides\t" + request.minUniquePeptides);
			pw.println("Input rows with a sequence\t" + r.inputRows);
			pw.println("Distinct sequences\t" + r.distinctSequences);
			pw.println("Sequences mapped to MAG proteins\t" + r.mappedSequences);
			pw.println("Genomes with mapped peptides\t" + r.genomesWithMappedPeptides);
			pw.println("Accepted genomes\t" + r.acceptedGenomes);
			pw.println("Sequences with razor protein\t" + r.peptidesWithRazorProtein);
			pw.println("Razor proteins\t" + r.analysis.proteins.length);
			pw.println("Genomes with intensity\t" + r.analysis.genomeIntensityMap.size());
		}
	}

	private static DecimalFormat intensityFormat() {
		return new DecimalFormat("0.###", DecimalFormatSymbols.getInstance(Locale.ROOT));
	}

	// ---------------------------------------------------------------- helpers

	private static Set<String> genomesOf(Iterable<String> proteins, String magIdentifier) {
		TreeSet<String> genomes = new TreeSet<>();
		for (String p : proteins) {
			if (isMag(p, magIdentifier))
				genomes.add(genomeOf(p));
		}
		return genomes;
	}

	private static boolean allMag(Iterable<String> proteins, String magIdentifier) {
		for (String p : proteins) {
			if (!isMag(p, magIdentifier))
				return false;
		}
		return true;
	}

	private static boolean isMag(String protein, String magIdentifier) {
		return protein != null && magIdentifier != null && !magIdentifier.isEmpty()
				&& protein.startsWith(magIdentifier);
	}

	/** Same genome rule as DdaMagPostAnalysisService / DdaMagProteinInferenceService. */
	private static String genomeOf(String protein) {
		int idx = protein.indexOf('_');
		return idx > 0 ? protein.substring(0, idx) : protein;
	}

	private static double sum(double[] values) {
		double t = 0;
		for (double v : values)
			t += v;
		return t;
	}

	private static MagDbItem resolveMagDb(String catalogueId, String version, File magdbConfigJson) {
		if (catalogueId == null || catalogueId.isBlank() || magdbConfigJson == null || !magdbConfigJson.isFile())
			return null;
		MagDbItem[][] all = MagDbConfigIO.parse(magdbConfigJson);
		if (all == null)
			return null;
		for (MagDbItem[] group : all) {
			if (group == null)
				continue;
			for (MagDbItem item : group) {
				if (item == null || !catalogueId.equals(item.getCatalogueID()))
					continue;
				if (version != null && !version.isBlank()) {
					for (MagDbItem.MagVersion mv : item.getVersionList()) {
						if (mv != null && version.equals(mv.getCatalogueVersion())) {
							item.setUsedVersion(mv);
							return item;
						}
					}
					return null; // an explicitly requested version that does not exist is an error
				}
				if (item.getUsedVersion() != null)
					return item;
			}
		}
		return null;
	}

	/** Peptide row of the input table; modSeq == sequence (no modification information is kept). */
	static final class TablePeptide extends MetaPeptide {
		TablePeptide(String sequence, String[] proteins, double[] intensity) {
			super(sequence, sequence, proteins, intensity, null);
			this.length = sequence.length();
		}

		@Override
		public Object[] getTableObjects() {
			return new Object[] { sequence, String.join(";", proteins), sum(intensity) };
		}

		@Override
		public Element getXmlPepElement() {
			throw new UnsupportedOperationException("TablePeptide is not serialised to XML");
		}
	}

	// ---------------------------------------------------------------- CLI

	private static final String USAGE = "Usage: PeptideTableMagHeadlessRunner --input <peptides.tsv|.csv> --out <dir>"
			+ " --catalogue <MAG catalogue id> [--version <catalogue version>] [--magdb-config <MAGDB_1_0.json>]"
			+ " [--func-db <func_def.db>] [--min-unique-peptides <N, default 0>] [--threads <n, default 4>]"
			+ " [--export-metaatlas]";

	public static void main(String[] args) throws Exception {
		Request r = new Request();
		for (int i = 0; i < args.length; i++) {
			String a = args[i];
			switch (a) {
			case "--input":
				r.peptideTable = new File(value(args, ++i, a));
				break;
			case "--out":
				r.outDir = new File(value(args, ++i, a));
				break;
			case "--catalogue":
				r.catalogueId = value(args, ++i, a);
				break;
			case "--version":
				r.catalogueVersion = value(args, ++i, a);
				break;
			case "--magdb-config":
				r.magdbConfigJson = new File(value(args, ++i, a));
				break;
			case "--func-db":
				r.funcDbFile = new File(value(args, ++i, a));
				break;
			case "--min-unique-peptides":
				r.minUniquePeptides = Integer.parseInt(value(args, ++i, a));
				break;
			case "--threads":
				r.threads = Integer.parseInt(value(args, ++i, a));
				break;
			case "--export-metaatlas":
				r.exportMetaAtlasParquet = true;
				break;
			case "-h":
			case "--help":
				System.out.println(USAGE);
				return;
			default:
				throw new IllegalArgumentException("Unknown argument: " + a + "\n" + USAGE);
			}
		}
		if (r.peptideTable == null || r.outDir == null || r.catalogueId == null)
			throw new IllegalArgumentException(USAGE);
		if (r.funcDbFile != null && !r.funcDbFile.isFile())
			throw new IllegalArgumentException("--func-db not found: " + r.funcDbFile);

		r.progressCallback = new ProgressCallback() {
			@Override
			public void onProgress(double workDone, double max) {
			}

			@Override
			public void onMessage(String message) {
				System.out.println(message);
			}

			@Override
			public void onComplete(Object result) {
				System.out.println("Complete: " + result);
			}

			@Override
			public void onError(Exception error) {
				error.printStackTrace(System.err);
			}
		};
		run(r);
	}

	private static String value(String[] args, int i, String flag) {
		if (i >= args.length)
			throw new IllegalArgumentException(flag + " requires a value");
		return args[i];
	}
}
