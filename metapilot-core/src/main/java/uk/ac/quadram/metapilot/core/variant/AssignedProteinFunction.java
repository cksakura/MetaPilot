package uk.ac.quadram.metapilot.core.variant;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashSet;
import java.util.Locale;
import java.util.Map;
import java.util.Set;
import java.util.TreeSet;
import java.util.function.Consumer;

import uk.ac.quadram.metapilot.task.mag.MagFuncSearcher;
import uk.ac.quadram.metapilot.task.mag.MetaProteinAnnoMag;

/**
 * Function of the searched-database proteins a substitution PSM was assigned
 * to, from the same per-genome eggNOG-mapper tables
 * ({@code <genome>_eggNOG.tsv}) the workflow's functional annotation reads.
 * <p>
 * A variant peptide is counted under those proteins everywhere else in the
 * workflow, and only species-representative genomes have eggNOG tables, so this
 * is the function the rest of the report already gives the PSM - copied beside
 * the substitution, not a new annotation of the variant's own carrier proteins.
 */
final class AssignedProteinFunction {

	private static final String SUFFIX = "_eggNOG.tsv";

	/** One protein's eggNOG row, the columns the candidates table carries. */
	static final class Annotation {
		final String keggKo;
		final String cog;
		final String cogCategory;
		final String preferredName;
		final String description;

		Annotation(String keggKo, String cog, String cogCategory, String preferredName, String description) {
			this.keggKo = keggKo;
			this.cog = cog;
			this.cogCategory = cogCategory;
			this.preferredName = preferredName;
			this.description = description;
		}
	}

	private final Map<String, Annotation> byProtein;
	/** Genomes whose eggNOG table was read; a protein of any other genome is NA, not unannotated. */
	private final Set<String> genomesRead;

	private AssignedProteinFunction(Map<String, Annotation> byProtein, Set<String> genomesRead) {
		this.byProtein = byProtein;
		this.genomesRead = genomesRead;
	}

	/**
	 * Reads the eggNOG rows of every protein the candidates were assigned to.
	 * Null when {@code eggnogDir} is null or not a directory.
	 */
	static AssignedProteinFunction load(File eggnogDir, Collection<VariantCandidate> candidates,
			Consumer<String> info, Consumer<String> warn) {
		if (eggnogDir == null || !eggnogDir.isDirectory()) {
			warn.accept("Pangenome variant check: no eggNOG folder" + (eggnogDir == null ? "" : " at " + eggnogDir)
					+ "; assigned_* function columns are NA");
			return null;
		}
		Map<String, Set<String>> proteinsByGenome = new HashMap<>();
		for (VariantCandidate c : candidates) {
			for (String protein : c.assignedProteins) {
				String id = accession(protein);
				String genome = genomeOf(id);
				if (genome != null) {
					proteinsByGenome.computeIfAbsent(genome, k -> new HashSet<>()).add(id);
				}
			}
		}
		Map<String, File> tables = new HashMap<>();
		File[] files = eggnogDir.listFiles();
		if (files != null) {
			for (File f : files) {
				String name = f.getName();
				if (f.isFile() && name.endsWith(SUFFIX)) {
					String genome = name.substring(0, name.length() - SUFFIX.length());
					int dot = genome.indexOf('.');
					tables.putIfAbsent(dot > 0 ? genome.substring(0, dot) : genome, f);
				}
			}
		}
		Map<String, Annotation> byProtein = new HashMap<>();
		Set<String> read = new HashSet<>();
		int missing = 0;
		for (Map.Entry<String, Set<String>> e : proteinsByGenome.entrySet()) {
			File table = tables.get(e.getKey());
			if (table == null) {
				missing++;
				continue;
			}
			try {
				readTable(table, e.getValue(), byProtein);
				read.add(e.getKey());
			} catch (IOException ex) {
				missing++;
				warn.accept("Pangenome variant check: cannot read " + table + ": " + ex.getMessage());
			}
		}
		info.accept(String.format(Locale.ROOT, "Pangenome variant check: assigned-protein function from %d eggNOG "
				+ "table(s) in %s, %,d proteins annotated%s", read.size(), eggnogDir, byProtein.size(),
				missing > 0 ? "; " + missing + " assigned genome(s) have no table, their proteins are NA" : ""));
		return new AssignedProteinFunction(byProtein, read);
	}

	private static void readTable(File table, Set<String> wanted, Map<String, Annotation> into) throws IOException {
		try (BufferedReader br = new BufferedReader(new FileReader(table, StandardCharsets.UTF_8))) {
			br.readLine();
			String line;
			while ((line = br.readLine()) != null) {
				if (line.startsWith("#")) {
					continue;
				}
				int tab = line.indexOf('\t');
				if (tab <= 0 || !wanted.contains(line.substring(0, tab))) {
					continue;
				}
				String[] f = line.split("\t", -1);
				into.put(f[0], new Annotation(field(f, MetaProteinAnnoMag.koID),
						MagFuncSearcher.cogAccession(field(f, MetaProteinAnnoMag.nogID)),
						field(f, MetaProteinAnnoMag.cogCatID), field(f, MetaProteinAnnoMag.prefnameID),
						field(f, MetaProteinAnnoMag.desID)));
			}
		}
	}

	/** Empty for a missing column or eggNOG's {@code -}. */
	private static String field(String[] f, int i) {
		String v = i < f.length ? f[i].trim() : "";
		return v.equals("-") ? "" : v;
	}

	/**
	 * Fills the candidate's {@code assigned*} function fields: the union over
	 * its assigned proteins. They stay null (NA) when none of those proteins'
	 * genomes has an eggNOG table; a protein whose table has no row for it is
	 * unannotated and adds nothing.
	 */
	void annotate(VariantCandidate c) {
		Set<String> ko = new TreeSet<>();
		Set<String> cog = new TreeSet<>();
		Set<String> category = new TreeSet<>();
		Set<String> name = new TreeSet<>();
		Set<String> description = new LinkedHashSet<>();
		boolean looked = false;
		for (String protein : c.assignedProteins) {
			String id = accession(protein);
			String genome = genomeOf(id);
			if (genome == null || !genomesRead.contains(genome)) {
				continue;
			}
			looked = true;
			Annotation a = byProtein.get(id);
			if (a == null) {
				continue;
			}
			split(a.keggKo, ko);
			split(a.cog, cog);
			if (!a.cogCategory.isEmpty()) {
				category.add(a.cogCategory);
			}
			if (!a.preferredName.isEmpty()) {
				name.add(a.preferredName);
			}
			if (!a.description.isEmpty()) {
				description.add(a.description);
			}
		}
		if (!looked) {
			return;
		}
		c.assignedKeggKo = String.join(";", ko);
		c.assignedCog = String.join(";", cog);
		c.assignedCogCategory = String.join(";", category);
		c.assignedPreferredName = String.join(";", name);
		c.assignedDescription = String.join(" | ", description);
	}

	private static void split(String values, Set<String> into) {
		for (String v : values.split(",")) {
			if (!v.isBlank()) {
				into.add(v.trim());
			}
		}
	}

	/** The accession before any description, as the FASTA header and the eggNOG table name it. */
	private static String accession(String protein) {
		String a = protein.trim();
		int space = a.indexOf(' ');
		return space > 0 ? a.substring(0, space) : a;
	}

	/** {@code MGYG000002463_01761} -&gt; {@code MGYG000002463}; null for a host or other non-MGYG accession. */
	private static String genomeOf(String accession) {
		int underscore = accession.lastIndexOf('_');
		if (underscore <= 0 || !accession.startsWith("MGYG")) {
			return null;
		}
		String genome = accession.substring(0, underscore);
		int dot = genome.indexOf('.');
		return dot > 0 ? genome.substring(0, dot) : genome;
	}
}
