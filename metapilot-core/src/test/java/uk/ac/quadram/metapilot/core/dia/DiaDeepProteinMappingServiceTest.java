package uk.ac.quadram.metapilot.core.dia;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.util.HashMap;
import java.util.HashSet;
import java.util.Set;

import org.junit.jupiter.api.Test;

class DiaDeepProteinMappingServiceTest {

	private static HashMap<String, HashSet<String>> map(String[]... entries) {
		HashMap<String, HashSet<String>> m = new HashMap<>();
		for (String[] e : entries) {
			HashSet<String> proteins = new HashSet<>();
			for (int i = 1; i < e.length; i++)
				proteins.add(e[i]);
			m.put(e[0], proteins);
		}
		return m;
	}

	@Test
	void genomesOfRoundEvidence() {
		HashMap<String, HashSet<String>> rounds = map(
				new String[] { "AAAK", "MGYG000000001_00001", "MGYG000000002_00010" },
				new String[] { "CCCR", "MGYG000000001_00002_42" },
				new String[] { "DDDK", "MGYG000000761.1_00003" });
		assertEquals(Set.of("MGYG000000001", "MGYG000000002", "MGYG000000761.1"),
				DiaDeepProteinMappingService.genomesOf(rounds));
	}

	@Test
	void markerEvidenceKeepsOnlySearchedGenomes() {
		HashMap<String, HashSet<String>> marker = map(
				new String[] { "SEARCHEDK", "MGYG000000001_00005" },
				new String[] { "MIXEDK", "MGYG000000001_00006", "MGYG000009999_00001" },
				new String[] { "NEVERSEARCHEDK", "MGYG000009999_00002", "MGYG000008888_00003" });
		DiaDeepProteinMappingService.restrictToGenomes(marker, Set.of("MGYG000000001", "MGYG000000002"));
		assertEquals(Set.of("MGYG000000001_00005"), marker.get("SEARCHEDK"));
		assertEquals(Set.of("MGYG000000001_00006"), marker.get("MIXEDK"));
		assertFalse(marker.containsKey("NEVERSEARCHEDK"));
		assertEquals(2, marker.size());
	}

	@Test
	void noSearchedGenomesKeepsNothing() {
		HashMap<String, HashSet<String>> marker = map(new String[] { "AAAK", "MGYG000000001_00001" });
		DiaDeepProteinMappingService.restrictToGenomes(marker, Set.of());
		assertTrue(marker.isEmpty());
	}
}
