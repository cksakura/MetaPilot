package uk.ac.quadram.metapilot.core.variant;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;

import java.util.ArrayList;
import java.util.List;

import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Test;

class VariantColumnDictionaryTest {

	@Test
	@DisplayName("Every column of every table is documented once, in order, and nothing else is")
	void dictionaryMatchesTheTables() {
		List<String> expected = new ArrayList<>();
		add(expected, PangenomeVariantCheck.CANDIDATES_FILE, PangenomeVariantCheck.CANDIDATE_COLUMNS);
		add(expected, PangenomeVariantCheck.SITES_FILE, VariantSiteTable.COLUMNS);
		add(expected, PangenomeVariantCheck.SUMMARY_FILE, VariantEvidenceSummary.COLUMNS);
		add(expected, PangenomeVariantCheck.SPECIES_FILE, VariantSpeciesSummary.COLUMNS);
		List<String> documented = new ArrayList<>();
		for (String[] e : VariantColumnDictionary.ENTRIES) {
			assertEquals(3, e.length);
			assertFalse(e[2].isBlank(), e[0] + " " + e[1]);
			if (!"*".equals(e[0])) {
				documented.add(e[0] + " " + e[1]);
			}
		}
		assertEquals(expected, documented);
	}

	private static void add(List<String> out, String file, String[] columns) {
		for (String c : columns) {
			out.add(file + " " + c);
		}
	}
}
