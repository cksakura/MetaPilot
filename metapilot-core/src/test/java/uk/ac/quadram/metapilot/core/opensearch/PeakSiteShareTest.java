package uk.ac.quadram.metapilot.core.opensearch;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.io.IOException;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Test;

class PeakSiteShareTest {

	private final PeakSiteShare.Config cfg = new PeakSiteShare.Config();

	@Test
	@DisplayName("A site is a substitution site only when residue plus shift is another amino acid")
	void substitutionSiteIsArithmetic() {
		assertTrue(PeakSiteShare.isSubstitutionSite('V', 14.01565, cfg), "V to L/I");
		assertTrue(PeakSiteShare.isSubstitutionSite('I', -14.01565, cfg), "I to V");
		assertTrue(PeakSiteShare.isSubstitutionSite('L', -14.01565, cfg), "L to V");
		assertTrue(PeakSiteShare.isSubstitutionSite('N', 0.98402, cfg), "N to D, the same mass as deamidation");
		assertFalse(PeakSiteShare.isSubstitutionSite('K', 14.01565, cfg), "methyl-K: no residue has that mass");
		assertFalse(PeakSiteShare.isSubstitutionSite('K', 42.010565, cfg), "acetyl-K");
		assertFalse(PeakSiteShare.isSubstitutionSite('K', 43.005814, cfg), "carbamyl-K");
		assertTrue(PeakSiteShare.isSubstitutionSite('A', 43.005814, cfg), "the same mass on A is A to N");
		assertFalse(PeakSiteShare.isSubstitutionSite('M', 15.994915, cfg), "oxidised M");
	}

	private static PeakSiteShare.PeakShare peak(PeakSiteShare.Result r, int id, int localised, int onSite) {
		PeakSiteShare.PeakShare p = new PeakSiteShare.PeakShare(id, id);
		p.localised = localised;
		p.onSubstitutionSite = onSite;
		r.peaks.put(id, p);
		return p;
	}

	@Test
	@DisplayName("A peak is withdrawn at the threshold share, and only with enough localised PSMs")
	void decision() {
		PeakSiteShare.Result r = new PeakSiteShare.Result();
		peak(r, 1, 100, 93);
		peak(r, 2, 100, 60);
		peak(r, 3, 100, 59);
		peak(r, 4, 19, 19);
		peak(r, 5, 20, 20);
		peak(r, 6, 0, 0);

		PeakSiteShare.decide(r, cfg);

		assertEquals(Set.of(1, 2, 5), r.withdrawnPeakIds());
		assertTrue(Double.isNaN(r.peaks.get(6).share()));
	}

	@Test
	@DisplayName("Sites outside the peptide do not count, and I and L count alike")
	void counting() {
		PeakSiteShare.PeakShare p = new PeakSiteShare.PeakShare(1, 14.0);
		PeakSiteShare.count(p, "AVK", 1, 14.01565, cfg);
		PeakSiteShare.count(p, "AVK", 2, 14.01565, cfg);
		PeakSiteShare.count(p, "AVK", 7, 14.01565, cfg);
		PeakSiteShare.count(p, "AVK", -1, 14.01565, cfg);

		assertEquals(2, p.localised);
		assertEquals(1, p.onSubstitutionSite);
	}

	@Test
	@DisplayName("Without spectra nothing is localised, so nothing is withdrawn, and the PSMs are reported")
	void noSpectra() throws IOException {
		OpenSearchPsm psm = new OpenSearchPsm();
		psm.setPeakId(3);
		psm.setHitRank(1);
		psm.setTarget(true);
		psm.setPeptide("ACDEFGHIK");
		psm.setFileName("sample");

		PeakSiteShare.Result r = PeakSiteShare.measure(List.of(psm), Set.of(3), new ShiftLocaliser(20), Map.of(), cfg,
				null);

		assertEquals(1, r.withoutSpectra);
		assertTrue(r.withdrawnPeakIds().isEmpty());
	}
}
