package uk.ac.quadram.metapilot.core.runner;

import static org.junit.jupiter.api.Assertions.assertEquals;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.junit.jupiter.api.Test;

import uk.ac.quadram.metapilot.core.mag.GenomeEvidenceTest.GenomeRow;
import uk.ac.quadram.metapilot.core.mag.GenomeEvidenceTest.PeptideClass;

class OpenSearchGenomeEvidenceTest {

	@Test
	void classifiesIsotopeErrorsAsUnmodified() {
		assertEquals(PeptideClass.UNMODIFIED, OpenSearchGenomeEvidence.classify(0.0, false));
		assertEquals(PeptideClass.UNMODIFIED, OpenSearchGenomeEvidence.classify(1.0034, false));
		assertEquals(PeptideClass.UNMODIFIED, OpenSearchGenomeEvidence.classify(2.0067 + 0.02, true));
		assertEquals(PeptideClass.UNEXPLAINED, OpenSearchGenomeEvidence.classify(3.0101, false));
		assertEquals(PeptideClass.UNEXPLAINED, OpenSearchGenomeEvidence.classify(0.05, false));
		assertEquals(PeptideClass.UNEXPLAINED, OpenSearchGenomeEvidence.classify(-1.0034, false));
		assertEquals(PeptideClass.KNOWN_SHIFT, OpenSearchGenomeEvidence.classify(15.9949, true));
		assertEquals(PeptideClass.KNOWN_SHIFT, OpenSearchGenomeEvidence.classify(Double.NaN, true));
	}

	private static GenomeRow row(String genome, boolean decoy, boolean accepted, int kept) {
		GenomeRow row = new GenomeRow();
		row.genome = genome;
		row.decoy = decoy;
		row.accepted = accepted;
		row.keptPeptides = kept;
		return row;
	}

	@Test
	void sameGenusAsNamesTheStrongerAcceptedGenomeOfTheGenus() {
		List<GenomeRow> rows = new ArrayList<>();
		rows.add(row("KA", false, true, 806));
		rows.add(row("KP", false, true, 538));
		rows.add(row("KV", false, true, 199));
		rows.add(row("KX", false, false, 3));
		rows.add(row("EC", false, true, 50));
		rows.add(row("NOLINEAGE", false, true, 900));
		rows.add(row("rev_KP", true, true, 4));
		Map<String, String> genus = new HashMap<>();
		genus.put("KA", "Klebsiella");
		genus.put("KP", "Klebsiella");
		genus.put("KV", "Klebsiella");
		genus.put("KX", "Klebsiella");
		genus.put("EC", "Escherichia");
		genus.put("rev_KP", "Klebsiella");

		Map<String, String> label = MagTaxonomyTsvExporter.sameGenusAs(rows, genus);
		Map<String, String> expected = new HashMap<>();
		expected.put("KP", "KA");
		expected.put("KV", "KA");
		assertEquals(expected, label);
	}
}
