package uk.ac.quadram.metapilot.dbSearch.pfind;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;

import java.util.HashMap;

import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Test;

import uk.ac.quadram.metapilot.task.dda.igc.par.DDAIgcPFindParameters;
import uk.ac.quadram.metapilot.task.dda.mag.par.DDAMagPfindParameters;

class PFindEnzymeTest {

	private static HashMap<String, PFindConfig.PFindItem> template() {
		HashMap<String, PFindConfig.PFindItem> map = new HashMap<>();
		map.put("[param]", new PFindConfig.PFindItem(1, true, "[param]", "", "", false));
		map.put("enzyme", new PFindConfig.PFindItem(2, false, "enzyme", "Trypsin KR _ C", "", false));
		map.put("digest", new PFindConfig.PFindItem(3, false, "digest", "3", "", false));
		return map;
	}

	@Test
	@DisplayName("pFind names keep pFind's meaning: Trypsin cleaves before proline, Trypsin_P does not")
	void trypsinPairUsesPFindSemantics() {
		assertEquals("Trypsin KR _ C", PFindEnzymes.definition("Trypsin"));
		assertEquals("Trypsin_P KR P C", PFindEnzymes.definition("trypsin_p"));
		assertEquals("Lys-K K _ C", PFindEnzymes.definition(" Lys-K "));
		assertEquals("Asp-N D _ N", PFindEnzymes.definition("Asp-N"));
	}

	@Test
	@DisplayName("Names pFind does not define are rejected instead of falling back to trypsin")
	void unknownNameIsRejected() {
		assertThrows(IllegalArgumentException.class, () -> PFindEnzymes.definition("Trypsin/P"));
		assertThrows(IllegalArgumentException.class, () -> PFindEnzymes.definition("GluC"));
		assertThrows(IllegalArgumentException.class, () -> PFindEnzymes.definition(null));
	}

	@Test
	@DisplayName("The selected enzyme replaces the template's enzyme line")
	void selectedEnzymeIsWritten() {
		HashMap<String, PFindConfig.PFindItem> map = template();
		PFindConfig.applyEnzyme(map, "Chymotrypsin");
		assertEquals("enzyme=Chymotrypsin FYWML _ C\n", map.get("enzyme").toString());
		assertEquals("digest=3\n", map.get("digest").toString());
	}

	@Test
	@DisplayName("A template without an enzyme line fails rather than searching an unchosen enzyme")
	void templateWithoutEnzymeLineFails() {
		HashMap<String, PFindConfig.PFindItem> map = template();
		map.remove("enzyme");
		assertThrows(IllegalStateException.class, () -> PFindConfig.applyEnzyme(map, "Trypsin"));
	}

	@Test
	@DisplayName("Defaults are the digest pFind searches always used, listed first for the dropdowns")
	void defaultsKeepThePreviousDigest() {
		assertEquals("Trypsin", PFindEnzymes.DEFAULT_ENZYME);
		assertEquals(PFindEnzymes.DEFAULT_ENZYME, PFindEnzymes.names().get(0));
		assertEquals(15, PFindEnzymes.names().size());
		assertEquals(PFindEnzymes.DEFAULT_ENZYME, new DDAMagPfindParameters().getEnzyme());
		assertEquals(PFindEnzymes.DEFAULT_ENZYME, new DDAIgcPFindParameters().getEnzyme());
	}
}
