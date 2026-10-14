package uk.ac.quadram.metapilot.core.variant;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Test;

class VariantEvidenceSummaryTest {

	private static VariantCandidate candidate(String base, int position, char to,
			VariantCandidate.Ms2Verdict verdict) {
		VariantCandidate c = new VariantCandidate();
		c.engine = VariantCandidate.Engine.OPEN_SEARCH;
		c.basePeptide = base;
		c.position = position;
		c.fromAa = base.charAt(position - 1);
		c.toAa = to;
		StringBuilder v = new StringBuilder(base);
		v.setCharAt(position - 1, to);
		c.variantPeptide = v.toString();
		c.ms2Verdict = verdict;
		return c;
	}

	/** Every peptide of the candidates, found in UHGP-100 when listed in {@code found}. */
	private static PangenomeVariantCheck.Evidence evidence(List<VariantCandidate> candidates, String... found) {
		Map<String, Integer> index = new LinkedHashMap<>();
		for (VariantCandidate c : candidates) {
			index.putIfAbsent(PeptideAutomaton.normalise(c.variantPeptide), index.size());
			index.putIfAbsent(PeptideAutomaton.normalise(c.basePeptide), index.size());
		}
		PangenomeVariantLookup.QueryHit[] hits = new PangenomeVariantLookup.QueryHit[index.size()];
		for (int i = 0; i < hits.length; i++) {
			hits[i] = new PangenomeVariantLookup.QueryHit();
		}
		for (String f : found) {
			hits[index.get(PeptideAutomaton.normalise(f))].proteins = 1;
		}
		return new PangenomeVariantCheck.Evidence(index, null, hits, null, PangenomeVariantSettings.SpeciesScope.ALL);
	}

	private static List<Map<String, String>> table(List<VariantEvidenceSummary.Row> rows) {
		List<Map<String, String>> out = new ArrayList<>();
		for (VariantEvidenceSummary.Row r : rows) {
			assertEquals(VariantEvidenceSummary.COLUMNS.length, r.values.length);
			Map<String, String> m = new HashMap<>();
			for (int i = 0; i < r.values.length; i++) {
				m.put(VariantEvidenceSummary.COLUMNS[i], r.values[i]);
			}
			out.add(m);
		}
		return out;
	}

	@Test
	@DisplayName("MS2 verdicts do not split the rows: groups, ALL, base")
	void verdictsDoNotSplitRows() {
		List<VariantCandidate> candidates = List.of(
				candidate("GASDK", 1, 'A', VariantCandidate.Ms2Verdict.SUBSTITUTION),
				candidate("GLSDK", 1, 'A', VariantCandidate.Ms2Verdict.MODIFICATION),
				candidate("PSTEK", 2, 'T', VariantCandidate.Ms2Verdict.SUBSTITUTION));
		List<Map<String, String>> rows = table(VariantEvidenceSummary.summarise(candidates, List.of(),
				evidence(candidates, "AASDK", "PTTEK", "GASDK", "GLSDK", "PSTEK"), true));

		assertFalse(Arrays.asList(VariantEvidenceSummary.COLUMNS).contains("ms2_verdict"));
		assertEquals(4, rows.size(), "two groups, one ALL, one base");
		assertEquals("G->A", rows.get(0).get("substitution"));
		assertEquals("2", rows.get(0).get("n_candidates"));
		assertEquals("50.00", rows.get(0).get("pct_in_pangenome_all"));
		assertEquals("ALL", rows.get(2).get("substitution"));
		assertEquals("3", rows.get(2).get("n_candidates"));
		assertEquals("base", rows.get(3).get("substitution"));
		assertEquals("substitution", rows.get(0).get("level"));
		assertEquals("1", rows.get(0).get("n_unique_in_pangenome_all"), "AASDK found, ALSDK not");
	}

	@Test
	@DisplayName("Open search: per-peak rows plus rows pooled over peaks")
	void peaksArePooledToo() {
		VariantCandidate a = candidate("GASDK", 1, 'A', VariantCandidate.Ms2Verdict.SUBSTITUTION);
		VariantCandidate b = candidate("GLSDK", 1, 'A', VariantCandidate.Ms2Verdict.SUBSTITUTION);
		a.peakId = 3;
		b.peakId = 7;
		a.peakSubstitutionShaped = true;
		b.peakSubstitutionShaped = true;
		List<VariantCandidate> candidates = List.of(a, b);
		List<Map<String, String>> rows = table(VariantEvidenceSummary.summarise(candidates, List.of(),
				evidence(candidates, "AASDK"), true));
		assertEquals(5, rows.size(), "pooled, two peaks, ALL, base");
		assertEquals("substitution", rows.get(0).get("level"));
		assertEquals("ALL", rows.get(0).get("peak_id"));
		assertEquals("2", rows.get(0).get("n_candidates"));
		assertEquals("substitution_peak", rows.get(1).get("level"));
		assertEquals("3", rows.get(1).get("peak_id"));
		assertEquals("7", rows.get(2).get("peak_id"));
		assertEquals("all", rows.get(3).get("level"));
		assertEquals("NA", rows.get(0).get("fisher_p_all"), "no controls, no test");
	}

	@Test
	@DisplayName("One-sided Fisher exact test and Benjamini-Hochberg")
	void fisherAndBh() {
		assertEquals(0.05, VariantEvidenceSummary.fisherGreater(3, 3, 0, 3), 1e-12, "1 / C(6,3)");
		assertEquals(1.0, VariantEvidenceSummary.fisherGreater(0, 3, 3, 3), 1e-12);
		// 2 of 4 found vs 1 of 4: P(X >= 2), X ~ hypergeometric(8, 3, 4) = (C(3,2)C(5,2) + C(3,3)C(5,1)) / C(8,4)
		assertEquals(35.0 / 70, VariantEvidenceSummary.fisherGreater(2, 4, 1, 4), 1e-12);
		assertTrue(Double.isNaN(VariantEvidenceSummary.fisherGreater(1, 1, 0, 0)));

		double[] q = VariantEvidenceSummary.benjaminiHochberg(new double[] { 0.01, Double.NaN, 0.04, 0.03 });
		assertEquals(0.03, q[0], 1e-12);
		assertTrue(Double.isNaN(q[1]));
		assertEquals(0.04, q[2], 1e-12);
		assertEquals(0.04, q[3], 1e-12);
	}
}
