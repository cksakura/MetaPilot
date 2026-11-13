package uk.ac.quadram.metapilot.core.dia;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertNull;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.io.File;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.util.List;

import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;

import uk.ac.quadram.metapilot.dbSearch.diann.DiannParameter;

class PeptideListLibrarySupportTest {

	@TempDir
	File dir;

	@Test
	void countsInternalTrypticSitesOnly() {
		assertEquals(0, PeptideListLibrarySupport.internalCleavageSites("AGLVLAEK"));
		assertEquals(1, PeptideListLibrarySupport.internalCleavageSites("AGKLVLAER"));
		assertEquals(2, PeptideListLibrarySupport.internalCleavageSites("AKGRLVLAEK"));
	}

	@Test
	void readsBaseSequencesStripsModificationsAndSorts() throws Exception {
		File table = new File(dir, "final_peptides.tsv");
		Files.write(table.toPath(), List.of(
				"Sequence\tBase sequence\tLength",
				"(Acetyl)MTEYK\tMTEYK\t5",
				"ZZPEPTIDEK\tPEPTIDEK\t8",
				"AAA(Oxidation)MK\tAAAMK\t5",
				"PEPTIDEK\tPEPTIDEK\t8"), StandardCharsets.UTF_8);
		PeptideListLibrarySupport.PeptideList list = PeptideListLibrarySupport.readPeptideTable(table);
		assertEquals(List.of("AAAMK", "MTEYK", "PEPTIDEK"), list.peptides);
		assertNull(PeptideListLibrarySupport.stripModifications("AA(Oxidation"));
		assertEquals("AAMK", PeptideListLibrarySupport.stripModifications("AA(Oxidation)MK"));
	}

	@Test
	void fastaIsReproducible() throws Exception {
		File table = new File(dir, "t.tsv");
		Files.write(table.toPath(), List.of("Base sequence", "TTTK", "AAAK", "GGGK"), StandardCharsets.UTF_8);
		PeptideListLibrarySupport.PeptideList list = PeptideListLibrarySupport.readPeptideTable(table);
		File f1 = new File(dir, "a.fasta");
		File f2 = new File(dir, "b.fasta");
		PeptideListLibrarySupport.writePeptideFasta(list, f1);
		PeptideListLibrarySupport.writePeptideFasta(PeptideListLibrarySupport.readPeptideFasta(f1), f2);
		assertEquals(Files.readAllLines(f1.toPath()), Files.readAllLines(f2.toPath()));
		assertEquals(List.of(">pro0 pro0", "AAAK", ">pro1 pro1", "GGGK", ">pro2 pro2", "TTTK"), Files.readAllLines(f1.toPath()));
	}

	@Test
	void peptideEntryRulesApplyOnlyDuringTheBuildStep() {
		DiannParameter p = new DiannParameter();
		p.setAcety_proN(true);
		p.setM_excision(true);
		boolean[] seen = new boolean[2];
		PeptideListLibrarySupport.runWithPeptideEntryRules(p, () -> {
			seen[0] = p.isAcety_proN();
			seen[1] = p.isM_excision();
		});
		assertFalse(seen[0]);
		assertFalse(seen[1]);
		assertTrue(p.isAcety_proN()); // restored for the later searches that reuse the parameter object
		assertTrue(p.isM_excision());
	}

	@Test
	void peptideLibraryRulesTurnOffProteinNTermAndKeepPeptidesIntact() throws Exception {
		File table = new File(dir, "t.tsv");
		Files.write(table.toPath(), List.of("Base sequence", "AGLVLAEK", "AKGRLVLAEK", "AKGRLKVLRAEK"), StandardCharsets.UTF_8);
		PeptideListLibrarySupport.PeptideList list = PeptideListLibrarySupport.readPeptideTable(table);
		DiannParameter p = new DiannParameter();
		p.setAcety_proN(true);
		p.setM_excision(true);
		p.setOxidation_M(true);
		int missed = PeptideListLibrarySupport.configureForPeptideLibrary(p, 1, list);
		assertEquals(2, missed); // most sites in the list is 4, capped at 2
		assertFalse(p.isAcety_proN());
		assertFalse(p.isM_excision());
		assertEquals(1, list.countWithMoreSitesThan(2));
		List<String> info = PeptideListLibrarySupport.describe(list, p);
		assertTrue(info.contains("Peptides with more than 2 internal cleavage sites (split):\t1"));
		String cmd = p.generateFasta(new File(dir, "x.fasta").getAbsolutePath(), new File(dir, "x.speclib").getAbsolutePath(), false);
		assertFalse(cmd.contains("UniMod:1,42.010565,*n"));
		assertFalse(cmd.contains("--met-excision"));
		assertTrue(cmd.contains("UniMod:35,15.994915,M"));
		assertTrue(cmd.contains("--missed-cleavages 2"));
	}
}
