package uk.ac.quadram.metapilot.core.runner;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;

import uk.ac.quadram.metapilot.core.opensearch.OpenSearchPsm;
import uk.ac.quadram.metapilot.core.opensearch.ShiftLocaliser;

class SetAsideReportTest {

	@TempDir
	File tmp;

	private static OpenSearchPsm psm(String file, int scan, String peptide, String protein, double delta,
			ShiftLocaliser.State state, String... explanations) {
		OpenSearchPsm p = new OpenSearchPsm();
		p.setFileName(file);
		p.setScan(scan);
		p.setPeptide(peptide);
		p.setTarget(!protein.startsWith("rev_"));
		p.setProtein(protein);
		p.setMassDiff(delta);
		p.setNormalizedMassDiff(delta);
		p.setPeakId(-1);
		p.setQValue(0.004);
		if (state != null) {
			p.setShiftLocalisation(new ShiftLocaliser.Result(5.0, 1, 0.9, state, Arrays.asList(explanations), delta));
		}
		return p;
	}

	private static List<String> rows(File f) throws IOException {
		List<String> out = new ArrayList<>();
		for (String line : Files.readAllLines(f.toPath())) {
			if (!line.startsWith("#")) {
				out.add(line);
			}
		}
		return out;
	}

	@Test
	@DisplayName("Reasons follow the localisation; peptides are grouped by reason and mass; genomes show what was "
			+ "set aside beside what was identified")
	void peptideAndGenomeTables() throws IOException {
		ShiftLocaliser.State explained = ShiftLocaliser.State.LOCALISED_EXPLAINED;
		OpenSearchPsm subst1 = psm("run1", 1, "SAMPLEK", "MGYG000000001_00010", 30.0106, explained, "subst A->T");
		OpenSearchPsm subst2 = psm("run2", 2, "SAMPLEK", "MGYG000000001_00010", 30.0101, explained, "subst A->T");
		OpenSearchPsm tie = psm("run1", 3, "PESTIDEK", "MGYG000000002_00020", 14.0157, explained, "subst S->T",
				"unimod Methyl");
		OpenSearchPsm unplaced = psm("run1", 4, "PEPTIDER", "MGYG000000002_00021", 4.0531,
				ShiftLocaliser.State.UNSUPPORTED);
		OpenSearchPsm decoy = psm("run1", 5, "KELPMAS", "rev_MGYG000000001_00010", 30.0106, explained, "subst A->T");
		OpenSearchPsm noSpectrum = psm("run1", 6, "PEPTIDEK", "MGYG000000002_00022", 7.7, null);
		assertEquals(SetAsideReport.Reason.SUBSTITUTION, SetAsideReport.reasonOf(subst1));
		assertEquals(SetAsideReport.Reason.MODIFICATION_OR_SUBSTITUTION, SetAsideReport.reasonOf(tie));
		assertEquals(SetAsideReport.Reason.SHIFT_NOT_ON_FRAGMENTS, SetAsideReport.reasonOf(unplaced));
		assertEquals(SetAsideReport.Reason.NOT_LOCALISED, SetAsideReport.reasonOf(noSpectrum));
		assertEquals("MGYG000000001", SetAsideReport.genomeOf("MGYG000000001_00010"));
		assertEquals("sp|P12345|HUMAN_X", SetAsideReport.genomeOf("sp|P12345|HUMAN_X"));

		List<OpenSearchPsm> identified = Arrays.asList(
				psm("run1", 10, "AAAAK", "MGYG000000001_00011", 0.0, null),
				psm("run1", 11, "AAAAK", "MGYG000000001_00011", 0.0, null),
				psm("run1", 12, "CCCCK", "MGYG000000001_00012", 0.0, null));
		int written = SetAsideReport.write(Arrays.asList(subst1, subst2, tie, unplaced, decoy, noSpectrum), identified,
				tmp, "rev_");
		assertEquals(4, written, "the two substitution PSMs of one peptide are one row; the decoy is left out");

		List<String> peptides = rows(new File(tmp, SetAsideReport.PEPTIDES_FILE));
		assertEquals(5, peptides.size());
		String[] first = peptides.get(1).split("\t", -1);
		assertEquals("SAMPLEK", first[0]);
		assertEquals("substitution", first[1]);
		assertEquals("2", first[3], "site, 1-based");
		assertEquals("subst A->T", first[4]);
		assertEquals("2", first[6], "PSMs");
		assertEquals("2", first[7], "samples");
		assertEquals("MGYG000000001", first[9]);

		List<String> genomes = rows(new File(tmp, SetAsideReport.GENOMES_FILE));
		assertEquals(3, genomes.size());
		String[] header = genomes.get(0).split("\t", -1);
		String[] g1 = genomes.get(1).split("\t", -1);
		assertEquals(header.length, g1.length);
		assertEquals("MGYG000000001", g1[0]);
		assertEquals("3", g1[1], "identified PSMs");
		assertEquals("2", g1[2], "identified peptides");
		assertEquals("2", g1[3], "set-aside PSMs");
		assertEquals("1", g1[4], "set-aside peptides");
		assertEquals("0.4000", g1[5]);
		String[] g2 = genomes.get(2).split("\t", -1);
		assertEquals("MGYG000000002", g2[0]);
		assertEquals("0", g2[1]);
		assertEquals("3", g2[3]);
		assertTrue(genomes.get(0).contains("modification_or_substitution_psms"));
	}
}
