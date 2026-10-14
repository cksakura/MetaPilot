package uk.ac.quadram.metapilot.core.variant;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertNotEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.io.File;
import java.io.IOException;
import java.io.OutputStream;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.zip.GZIPOutputStream;

import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;

class PangenomeVariantCheckTest {

	@TempDir
	File tmp;

	private static VariantCandidate candidate(String file, int scan, String base, int position, char to) {
		VariantCandidate c = new VariantCandidate();
		c.engine = VariantCandidate.Engine.PFIND;
		c.file = file;
		c.scan = scan;
		c.psmId = 2;
		c.basePeptide = base;
		c.position = position;
		c.fromAa = base.charAt(position - 1);
		c.toAa = to;
		StringBuilder v = new StringBuilder(base);
		v.setCharAt(position - 1, to);
		c.variantPeptide = v.toString();
		return c;
	}

	private void write(File f, String content, boolean gzip) throws IOException {
		f.getParentFile().mkdirs();
		if (gzip) {
			try (OutputStream out = new GZIPOutputStream(Files.newOutputStream(f.toPath()))) {
				out.write(content.getBytes(StandardCharsets.US_ASCII));
			}
		} else {
			Files.write(f.toPath(), content.getBytes(StandardCharsets.US_ASCII));
		}
	}

	/** A miniature MAG version directory with a prepared UHGP-100 resource. */
	private File magVersionDir() throws IOException {
		File version = new File(tmp, "human-gut/v2.0.2");
		write(new File(version, PangenomeResource.METADATA), "Genome\tSpecies_rep\n"
				+ "MGYG000000001\tMGYG000000001\n" + "MGYG000000002\tMGYG000000002\n"
				+ "MGYG000000003\tMGYG000000003\n", false);
		File catalogue = new File(version, PangenomeResource.DEFAULT_SUBDIR);
		write(new File(catalogue, PangenomeResource.FAA),
				">MGYG000000001_00001\nKPEPTIDEKAAAVEKR\n" // base PEPTIDEK, variant PEPTLDEK? no
						+ ">MGYG000000002_00001\nKPEPSIDEKGGG\n" // PEPSIDEK: variant of PEPTIDEK, species 2
						+ ">MGYG000000003_00001\nRAAAGEKR\n", // AAAGEK: variant of AAAVEK, species 3
				true);
		write(new File(catalogue, PangenomeResource.REP_SPECIES), "", true);
		return version;
	}

	private PangenomeVariantCheck.Request request(File version, List<VariantCandidate> candidates) throws IOException {
		File searched = new File(tmp, "searched.fasta");
		write(searched, ">sp1\nKPEPTIDEK\n>rev_sp2\nKPEPSIDEK\n>sp3\nKAAAVEK\n", false);
		File taxonomy = new File(tmp, "taxonomy_analysis");
		write(new File(taxonomy, "genome_final_selected.tsv"),
				"Genome\tReason\nMGYG000000002\tGREEDY_core\n\n#PrivilegeGenomes\t0\n", false);
		PangenomeVariantCheck.Request r = new PangenomeVariantCheck.Request();
		r.engine = VariantCandidate.Engine.PFIND;
		r.candidates = candidates;
		r.identifiedBasePeptides = List.of("PEPTIDEK", "AAAVEK", "GGSTTK");
		r.searchedFastas = List.of(searched);
		r.decoyPrefix = "rev_";
		r.magVersionDir = version;
		r.taxonomyDir = taxonomy;
		r.outputDir = new File(tmp, "out");
		r.threads = 2;
		r.undecidableToleranceDa = 0.005;
		return r;
	}

	private static List<Map<String, String>> readTsv(File f) throws IOException {
		List<String> lines = Files.readAllLines(f.toPath(), StandardCharsets.UTF_8);
		String[] header = lines.get(0).split("\t", -1);
		List<Map<String, String>> rows = new ArrayList<>();
		for (String line : lines.subList(1, lines.size())) {
			String[] v = line.split("\t", -1);
			assertEquals(header.length, v.length, line);
			Map<String, String> row = new HashMap<>();
			for (int i = 0; i < header.length; i++) {
				row.put(header[i], v[i]);
			}
			rows.add(row);
		}
		return rows;
	}

	@Test
	@DisplayName("End to end: classes, scopes and summary on a miniature catalogue")
	void endToEnd() throws IOException {
		File version = magVersionDir();
		List<VariantCandidate> candidates = new ArrayList<>();
		VariantCandidate inPangenome = candidate("raw1", 10, "PEPTIDEK", 4, 'S');
		VariantCandidate novel = candidate("raw1", 11, "PEPTIDEK", 2, 'W');
		VariantCandidate otherSpecies = candidate("raw2", 5, "AAAVEK", 4, 'G');
		VariantCandidate undecidable = candidate("raw2", 6, "AAAVEK", 4, 'G');
		undecidable.ms2Decidable = false;
		undecidable.isobaricModifications.add("test alias");
		inPangenome.assignedProteins.add("MGYG000000002_00007");
		otherSpecies.assignedProteins.add("MGYG000000001_00003");
		otherSpecies.assignedProteins.add("HUMAN_P12345");
		candidates.add(novel);
		candidates.add(inPangenome);
		candidates.add(otherSpecies);
		candidates.add(undecidable);

		PangenomeVariantCheck.Outcome outcome = PangenomeVariantCheck.run(request(version, candidates));
		assertTrue(outcome.ran, outcome.message);

		assertEquals(VariantCandidate.PangenomeClass.SUBSTITUTION_IN_PANGENOME, inPangenome.pangenomeClass);
		assertEquals(Boolean.FALSE, inPangenome.inSearchedDb, "the searched DB has it only as a decoy");
		assertEquals(Boolean.TRUE, inPangenome.inPangenomeDetected);
		assertEquals(1, inPangenome.nSpecies);
		assertEquals("MGYG000000002_00001", inPangenome.exampleProtein);

		assertEquals(VariantCandidate.PangenomeClass.SUBSTITUTION_NOVEL, novel.pangenomeClass);
		assertFalse(novel.inPangenomeAll);

		assertEquals(VariantCandidate.PangenomeClass.SUBSTITUTION_IN_PANGENOME, otherSpecies.pangenomeClass);
		assertEquals(Boolean.FALSE, otherSpecies.inPangenomeDetected, "species 3 was not selected");

		assertEquals(Boolean.TRUE, inPangenome.inSameSpecies, "variant found in species 2, protein from species 2");
		assertEquals(Boolean.FALSE, otherSpecies.inSameSpecies, "variant only in species 3, protein from species 1");
		assertEquals(java.util.Set.of("MGYG000000001"), otherSpecies.assignedSpecies, "host accession ignored");
		assertEquals(null, novel.inSameSpecies, "no assigned protein");

		assertEquals(VariantCandidate.PangenomeClass.UNDECIDABLE_BY_MS2, undecidable.pangenomeClass);
		assertTrue(undecidable.inPangenomeAll, "pangenome columns still carry the evidence");

		List<Map<String, String>> rows = readTsv(outcome.candidatesFile);
		assertEquals(4, rows.size());
		assertEquals("raw1", rows.get(0).get("file"));
		assertEquals("10", rows.get(0).get("scan"), "rows are sorted by file, scan, psm_id, peptide, position");
		assertEquals("NA", rows.get(0).get("hit_rank"), "pFind has no hit rank");
		assertEquals("MGYG000000002_00007", rows.get(0).get("assigned_proteins"));
		assertEquals("NA", rows.get(0).get("protein_sites"), "not in the searched database");
		assertEquals("1", rows.get(0).get("n_substitutions_in_variant_peptide"));
		assertEquals("1", rows.get(0).get("n_candidates_for_psm"));
		assertEquals("NA", rows.get(0).get("localized_sequence"));
		assertEquals("NA", rows.get(0).get("reassigner_in_db"));
		assertTrue(outcome.sitesFile.isFile());
		assertTrue(outcome.columnsFile.isFile());
		assertEquals("NA", rows.get(0).get("peak_id"));
		assertEquals("pfind", rows.get(0).get("engine"));
		assertEquals("NOT_APPLICABLE", rows.get(0).get("ms2_verdict"), "pFind calls are not re-ranked");
		assertEquals("NA", rows.get(0).get("ms2_score_substitution"));
		assertEquals("NA", rows.get(0).get("ms2_score_elsewhere"));
		assertEquals("NA", rows.get(0).get("ms2_localises_elsewhere"));
		assertEquals("NA", rows.get(0).get("ms2_best_form"));
		assertEquals("true", rows.get(0).get("in_same_species"), "raw1 scan 10 is the in-pangenome candidate");
		assertEquals("MGYG000000002", rows.get(0).get("assigned_species"));
		assertEquals("NA", rows.get(1).get("in_same_species"));
		assertEquals("NA", rows.get(0).get("n_genomes"), "no membership table, no genome counts");
		assertEquals("NA", rows.get(0).get("example_genome_same_species"));

		List<Map<String, String>> summary = readTsv(outcome.summaryFile);
		Map<String, String> base = summary.get(summary.size() - 1);
		assertEquals("base", base.get("substitution"));
		assertEquals("100.00", base.get("pct_in_pangenome_all"));
		boolean sawAll = false;
		for (Map<String, String> row : summary) {
			if ("ALL".equals(row.get("substitution")) && "true".equals(row.get("ms2_decidable"))) {
				sawAll = true;
				assertEquals("3", row.get("n_candidates"));
				assertEquals("66.67", row.get("pct_in_pangenome_all"));
				assertEquals("33.33", row.get("pct_in_pangenome_detected"));
			}
		}
		assertTrue(sawAll);
	}

	/** Genome types and a cluster membership on top of {@link #magVersionDir()}. */
	private File magVersionDirWithMembership() throws IOException {
		File version = magVersionDir();
		write(new File(version, PangenomeResource.METADATA), "Genome\tGenome_type\tSpecies_rep\n"
				+ "MGYG000000001\tIsolate\tMGYG000000001\n" + "MGYG000000002\tIsolate\tMGYG000000002\n"
				+ "MGYG000000003\tMAG\tMGYG000000003\n" + "MGYG000000004\tMAG\tMGYG000000002\n"
				+ "MGYG000000005\tMAG\tMGYG000000002\n" + "MGYG000000006\tIsolate\tMGYG000000003\n", false);
		write(new File(version, PangenomeResource.DEFAULT_SUBDIR + "/uhgp-100/" + UhgpMembership.FILE),
				"MGYG000000001_00001\tMGYG000000001_00001\n" + "MGYG000000002_00001\tMGYG000000002_00001\n"
						+ "MGYG000000002_00001\tMGYG000000004_00010\n"
						+ "MGYG000000002_00001\tMGYG000000006_00003\n"
						+ "MGYG000000003_00001\tMGYG000000003_00001",
				false);
		return version;
	}

	@Test
	@DisplayName("Species summary and assigned-protein function")
	void speciesSummaryAndAssignedFunction() throws IOException {
		File version = magVersionDirWithMembership();
		File eggnog = new File(version, "eggNOG");
		write(new File(eggnog, "MGYG000000002_eggNOG.tsv"),
				"#query\tseed_ortholog\tevalue\tscore\teggNOG_OGs\tmax_annot_lvl\tCOG_category\tDescription"
						+ "\tPreferred_name\tGOs\tEC\tKEGG_ko\n"
						+ "MGYG000000002_00007\tx\t0\t1\tCOG0001@1|root,COG0002@2|Bacteria,4NE6@976|Bacteroidetes"
						+ "\t2|Bacteria\tE\tAmino acid permease\taroA\t-\t-\tko:K00800,ko:K00801\n"
						+ "MGYG000000002_00099\tx\t0\t1\t-\t-\tJ\tunrelated\t-\t-\t-\t-\n",
				false);
		VariantCandidate inPangenome = candidate("raw1", 10, "PEPTIDEK", 4, 'S');
		inPangenome.assignedProteins.add("MGYG000000002_00007");
		inPangenome.assignedProteins.add("MGYG000000002_00008");
		VariantCandidate otherSpecies = candidate("raw2", 5, "AAAVEK", 4, 'G');
		otherSpecies.assignedProteins.add("MGYG000000001_00003");
		VariantCandidate unannotated = candidate("raw2", 7, "AAAVEK", 4, 'G');
		unannotated.assignedProteins.add("MGYG000000002_00008");

		PangenomeVariantCheck.Request request = request(version,
				new ArrayList<>(List.of(inPangenome, otherSpecies, unannotated)));
		request.eggnogDir = eggnog;
		PangenomeVariantCheck.Outcome outcome = PangenomeVariantCheck.run(request);
		assertTrue(outcome.ran, outcome.message);

		List<Map<String, String>> rows = readTsv(outcome.candidatesFile);
		Map<String, String> annotated = rows.get(0);
		assertEquals("ko:K00800;ko:K00801", annotated.get("assigned_kegg_ko"));
		assertEquals("COG0001;COG0002", annotated.get("assigned_cog"), "root and Bacteria COG, as the report has");
		assertEquals("E", annotated.get("assigned_cog_category"));
		assertEquals("aroA", annotated.get("assigned_preferred_name"));
		assertEquals("Amino acid permease", annotated.get("assigned_description"), "_00008 has no row, adds nothing");
		assertEquals("NA", rows.get(1).get("assigned_kegg_ko"), "species 1 has no eggNOG table");
		assertEquals("", rows.get(2).get("assigned_kegg_ko"), "table read, protein not annotated");
		assertEquals("", rows.get(2).get("assigned_description"));

		List<Map<String, String>> species = readTsv(outcome.speciesFile);
		assertEquals(2, species.size());
		Map<String, String> s2 = species.get(0);
		assertEquals("MGYG000000002", s2.get("species"), "selected species first");
		assertEquals("true", s2.get("in_genome_final_selected"));
		assertEquals("3", s2.get("n_genomes_in_species"));
		assertEquals("1", s2.get("n_isolates_in_species"));
		assertEquals("2", s2.get("n_candidates"));
		assertEquals("0", s2.get("n_candidates_shared"));
		assertEquals("2", s2.get("n_candidates_in_pangenome"));
		assertEquals("0", s2.get("n_candidates_novel"));
		assertEquals("2", s2.get("n_peptides"));
		assertEquals("2", s2.get("n_peptides_in_pangenome"));
		assertEquals("1", s2.get("n_peptides_in_pangenome_this_species"), "PEPSIDEK");
		assertEquals("1", s2.get("n_peptides_in_pangenome_other_species_only"), "AAAGEK, species 3 only");
		assertEquals("1", s2.get("n_other_genomes_carrying"), "genome 4; the representative 2 and species-3 genome 6 "
				+ "are not counted");
		assertEquals("0", s2.get("n_other_isolates_carrying"));

		Map<String, String> s1 = species.get(1);
		assertEquals("MGYG000000001", s1.get("species"));
		assertEquals("false", s1.get("in_genome_final_selected"));
		assertEquals("1", s1.get("n_genomes_in_species"));
		assertEquals("1", s1.get("n_candidates"));
		assertEquals("1", s1.get("n_peptides_in_pangenome"));
		assertEquals("0", s1.get("n_peptides_in_pangenome_this_species"));
		assertEquals("1", s1.get("n_peptides_in_pangenome_other_species_only"));
		assertEquals("0", s1.get("n_other_genomes_carrying"));
	}

	@Test
	@DisplayName("A PSM shared by two species counts under both, flagged as shared")
	void sharedCandidateCountsUnderEverySpecies() throws IOException {
		File version = magVersionDirWithMembership();
		VariantCandidate shared = candidate("raw1", 10, "PEPTIDEK", 4, 'S');
		shared.assignedProteins.add("MGYG000000001_00002");
		shared.assignedProteins.add("MGYG000000002_00007");
		PangenomeVariantCheck.Outcome outcome = PangenomeVariantCheck
				.run(request(version, new ArrayList<>(List.of(shared))));
		assertTrue(outcome.ran, outcome.message);

		List<Map<String, String>> species = readTsv(outcome.speciesFile);
		assertEquals(2, species.size());
		for (Map<String, String> row : species) {
			assertEquals("1", row.get("n_candidates"), row.get("species"));
			assertEquals("1", row.get("n_candidates_shared"), row.get("species"));
		}
		assertEquals("1", species.get(0).get("n_peptides_in_pangenome_this_species"), "PEPSIDEK is species 2");
		assertEquals("0", species.get(1).get("n_peptides_in_pangenome_this_species"));
		assertEquals("NA", readTsv(outcome.candidatesFile).get(0).get("assigned_kegg_ko"), "no eggNOG folder given");
	}

	@Test
	@DisplayName("Genome counts from the cluster membership: same species, species total, isolates, example")
	void genomeCounts() throws IOException {
		File version = magVersionDirWithMembership();
		VariantCandidate inPangenome = candidate("raw1", 10, "PEPTIDEK", 4, 'S');
		inPangenome.assignedProteins.add("MGYG000000002_00007");
		VariantCandidate otherSpecies = candidate("raw2", 5, "AAAVEK", 4, 'G');
		otherSpecies.assignedProteins.add("MGYG000000001_00003");
		VariantCandidate novel = candidate("raw1", 11, "PEPTIDEK", 2, 'W');

		PangenomeVariantCheck.Outcome outcome = PangenomeVariantCheck
				.run(request(version, new ArrayList<>(List.of(inPangenome, otherSpecies, novel))));
		assertTrue(outcome.ran, outcome.message);

		assertEquals(3, inPangenome.nGenomes, "genomes 2, 4 and 6 carry the PEPSIDEK cluster");
		assertEquals(2, inPangenome.nGenomesSameSpecies, "2 and 4 are species 2; 6 is species 3");
		assertEquals(3, inPangenome.nGenomesInSpecies, "species 2 has genomes 2, 4 and 5");
		assertEquals(1, inPangenome.nIsolatesSameSpecies);
		assertEquals("MGYG000000002", inPangenome.exampleGenomeSameSpecies);

		assertEquals(1, otherSpecies.nGenomes);
		assertEquals(0, otherSpecies.nGenomesSameSpecies);
		assertEquals(1, otherSpecies.nGenomesInSpecies);
		assertEquals("", otherSpecies.exampleGenomeSameSpecies);

		assertEquals(0, novel.nGenomes);
		assertEquals(-1, novel.nGenomesSameSpecies, "no assigned protein, species unknown");

		List<Map<String, String>> rows = readTsv(outcome.candidatesFile);
		assertEquals("3", rows.get(0).get("n_genomes"));
		assertEquals("2", rows.get(0).get("n_genomes_same_species"));
		assertEquals("3", rows.get(0).get("n_genomes_in_species"));
		assertEquals("1", rows.get(0).get("n_isolates_same_species"));
		assertEquals("MGYG000000002", rows.get(0).get("example_genome_same_species"));
		assertEquals("NA", rows.get(1).get("n_genomes_same_species"));
	}

	@Test
	@DisplayName("Protein sites: located in the searched database, rolled up per site, multi-site PSMs counted")
	void proteinSites() throws IOException {
		File version = magVersionDir();
		VariantCandidate a = candidate("raw1", 10, "PEPTIDEK", 4, 'S');
		a.assignedProteins.add("sp1");
		a.assignedProteins.add("sp9");
		a.qValue = 0.004;
		VariantCandidate b = candidate("raw2", 3, "PEPTIDEK", 4, 'S');
		b.assignedProteins.add("sp1");
		b.qValue = 0.001;
		VariantCandidate first = candidate("raw1", 12, "AAAVEK", 1, 'G');
		VariantCandidate second = candidate("raw1", 12, "AAAVEK", 4, 'G');
		first.variantPeptide = "GAAGEK";
		second.variantPeptide = "GAAGEK";
		first.assignedProteins.add("sp3");
		second.assignedProteins.add("sp3");

		PangenomeVariantCheck.Outcome outcome = PangenomeVariantCheck
				.run(request(version, new ArrayList<>(List.of(a, b, first, second))));
		assertTrue(outcome.ran, outcome.message);

		Map<String, Map<String, String>> byKey = new HashMap<>();
		for (Map<String, String> row : readTsv(outcome.candidatesFile)) {
			byKey.put(row.get("file") + ":" + row.get("scan") + ":" + row.get("position"), row);
		}
		Map<String, String> ra = byKey.get("raw1:10:4");
		assertEquals("sp1;sp9", ra.get("assigned_proteins"));
		assertEquals("sp1:5", ra.get("protein_sites"), "KPEPTIDEK: T is residue 5; sp9 is not in the database");
		assertEquals("K", ra.get("residue_before"));
		assertEquals("-", ra.get("residue_after"));
		assertEquals("4.000e-03", ra.get("psm_q_value"));
		Map<String, String> r1 = byKey.get("raw1:12:1");
		assertEquals("sp3:2", r1.get("protein_sites"));
		assertEquals("2", r1.get("n_substitutions_in_variant_peptide"));
		assertEquals("2", r1.get("n_candidates_for_psm"));
		assertEquals("sp3:5", byKey.get("raw1:12:4").get("protein_sites"));

		List<Map<String, String>> sites = readTsv(outcome.sitesFile);
		assertEquals(4, sites.size());
		Map<String, String> sp1 = sites.get(0);
		assertEquals("sp1", sp1.get("protein"));
		assertEquals("5", sp1.get("protein_position"));
		assertEquals("T->S", sp1.get("substitution"));
		assertEquals("true", sp1.get("located"));
		assertEquals("2", sp1.get("n_candidates"));
		assertEquals("1", sp1.get("n_candidates_shared"));
		assertEquals("2", sp1.get("n_psms"));
		assertEquals("raw1;raw2", sp1.get("files"));
		assertEquals("1.000e-03", sp1.get("min_psm_q_value"));
		assertEquals("PEPTIDEK@4", sp1.get("peptide_sites"));
		assertEquals("PEPSIDEK", sp1.get("variant_peptides"));
		assertEquals("true", sp1.get("in_pangenome_all"));
		assertEquals("NA", sp1.get("n_ms2_ranked"), "pFind is not ranked");
		assertEquals("2", sites.get(1).get("protein_position"));
		assertEquals("sp3", sites.get(1).get("protein"));
		assertEquals("5", sites.get(2).get("protein_position"));
		Map<String, String> sp9 = sites.get(3);
		assertEquals("sp9", sp9.get("protein"));
		assertEquals("NA", sp9.get("protein_position"));
		assertEquals("false", sp9.get("located"));
		assertEquals("1", sp9.get("n_candidates"));
	}

	@Test
	@DisplayName("Off, or no resource: nothing is written")
	void offOrAbsent() throws IOException {
		File version = magVersionDir();
		PangenomeVariantCheck.Request off = request(version, List.of(candidate("raw1", 1, "PEPTIDEK", 4, 'S')));
		off.settings.setMode(PangenomeVariantSettings.Mode.OFF);
		assertFalse(PangenomeVariantCheck.run(off).ran);
		assertFalse(PangenomeVariantCheck.wouldRun(off.settings, version));

		PangenomeVariantCheck.Request absent = request(new File(tmp, "no-such-version"),
				List.of(candidate("raw1", 1, "PEPTIDEK", 4, 'S')));
		absent.settings.setMode(PangenomeVariantSettings.Mode.ON);
		assertFalse(PangenomeVariantCheck.run(absent).ran);

		assertFalse(new File(tmp, "out").exists());
	}

	@Test
	@DisplayName("Controls are reproducible for a seed and change with it")
	void controlsDeterministic() {
		List<VariantCandidate> candidates = new ArrayList<>();
		for (int i = 0; i < 50; i++) {
			candidates.add(candidate("raw", i, "PEPTIDEK", 4, 'S'));
			candidates.add(candidate("raw", i, "AALGEK", 3, 'V'));
		}
		List<String> pool = List.of("GGTTK", "TAAK", "LLIIK", "MTLK", "PPPK", "ITTTR");
		List<String> a = peptides(MatchedControlGenerator.generate(candidates, pool, 1, 2));
		List<String> b = peptides(MatchedControlGenerator.generate(candidates, new ArrayList<>(pool), 1, 2));
		List<String> c = peptides(MatchedControlGenerator.generate(candidates, pool, 2, 2));
		assertEquals(200, a.size());
		assertEquals(a, b);
		assertNotEquals(a, c);
		for (MatchedControlGenerator.Control k : MatchedControlGenerator.generate(candidates, pool, 1, 1)) {
			char from = Character.toUpperCase(k.basePeptide.charAt(k.position - 1));
			char expected = k.parent.fromAa == 'I' ? 'L' : k.parent.fromAa;
			assertEquals(expected, from == 'I' ? 'L' : from, "control site carries the candidate's from-residue");
			assertEquals(k.parent.toAa, k.peptide.charAt(k.position - 1));
		}
	}

	private static List<String> peptides(List<MatchedControlGenerator.Control> controls) {
		List<String> out = new ArrayList<>();
		for (MatchedControlGenerator.Control c : controls) {
			out.add(c.peptide);
		}
		return out;
	}
}
