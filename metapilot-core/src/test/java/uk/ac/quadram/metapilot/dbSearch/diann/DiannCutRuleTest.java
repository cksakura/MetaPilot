package uk.ac.quadram.metapilot.dbSearch.diann;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Test;

class DiannCutRuleTest {

	private static DiannParameter withEnzyme(String enzyme) {
		DiannParameter par = new DiannParameter();
		par.setEnzyme(enzyme);
		return par;
	}

	@Test
	@DisplayName("Enzyme names from the DIA panels resolve regardless of case")
	void enzymeNamesFromThePanelsResolve() {
		assertEquals("K*,R*", withEnzyme("Trypsin/P").getCut());
		assertEquals("K*,R*,!*P", withEnzyme("Trypsin").getCut());
		assertEquals("K*", withEnzyme("Lys-C").getCut());
		assertEquals("F*,Y*,W*,M*,L*,!*P", withEnzyme("Chymotrypsin").getCut());
		assertEquals("E*", withEnzyme("GluC").getCut());
	}

	@Test
	@DisplayName("Asp-N cleaves before D")
	void aspNCutsBeforeAspartate() {
		assertEquals("*D", withEnzyme("AspN").getCut());
	}

	@Test
	@DisplayName("The DIA-NN command line carries the cut rule, and an unknown enzyme stops the run")
	void commandLineCutRule() {
		assertTrue(withEnzyme("Chymotrypsin").generateFasta("in.fasta", "lib.speclib", false)
				.contains("--cut F*,Y*,W*,M*,L*,!*P "));
		assertTrue(withEnzyme("AspN").generateFasta("in.fasta", "lib.speclib", false).contains("--cut *D "));
		assertThrows(IllegalArgumentException.class,
				() -> withEnzyme("Pepsin").generateFasta("in.fasta", "lib.speclib", false));
	}
}
