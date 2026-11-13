package uk.ac.quadram.metapilot.task.mag;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNull;

import java.io.File;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.util.List;

import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;

import uk.ac.quadram.metapilot.dbSearch.diann.DiaNNPrecursor;

class MarkerPeptideRatioTest {

	@TempDir
	File dir;

	private static DiaNNPrecursor precursor(String seq, String... proteins) {
		return new DiaNNPrecursor("run", proteins, seq, seq, 2, 0, 0, 0.001, 0, 0.001, 0.001, 0.001, 0, 0);
	}

	@Test
	void nonMarkerProteinsOfAMarkerGenomeDoNotCount() throws Exception {
		File markerFasta = new File(dir, "Marker4.fasta");
		Files.write(markerFasta.toPath(), List.of(
				">MGYG000000001_00010 DNA polymerase III subunit alpha", "PEPTIDEK",
				">MGYG000000001_00068 Valine--tRNA ligase", "PEPTIDER"), StandardCharsets.UTF_8);
		List<DiaNNPrecursor> precursors = List.of(
				precursor("AAAAK", "MGYG000000001_00010"),           // marker protein
				precursor("CCCCK", "MGYG000000001_00068_1"),         // marker protein with a suffix
				precursor("DDDDK", "MGYG000000001_00001"),           // same genome, not a marker
				precursor("EEEEK", "MGYG000000001_00002"),           // same genome, not a marker
				precursor("FFFFK", "MGYG000000002_00010"),           // other genome
				precursor("GGGGK", "REV_MGYG000000001_00010"));      // decoy
		MarkerGenomeSelector.MarkerPeptideRatioResult r =
				MarkerGenomeSelector.computeMarkerPeptideRatio(precursors, markerFasta);
		assertEquals(6, r.totalUniquePeptides);
		assertEquals(2, r.markerUniquePeptides); // before the fix: 4 (every MGYG000000001 protein matched)
		assertEquals(2.0 / 6, r.ratio, 1e-9);
	}

	@Test
	void suffixIsStrippedOnlyWhenAProteinIdRemains() {
		assertEquals("MGYG000000001_00010", MarkerGenomeSelector.stripLastSuffixKeepingProteinId("MGYG000000001_00010_1"));
		assertNull(MarkerGenomeSelector.stripLastSuffixKeepingProteinId("MGYG000000001_00010"));
		assertNull(MarkerGenomeSelector.stripLastSuffixKeepingProteinId("MGYG000000001"));
	}
}
