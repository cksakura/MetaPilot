package uk.ac.quadram.metapilot.core.opensearch;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertNotNull;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.nio.charset.StandardCharsets;
import java.util.List;

import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Test;

import uk.ac.quadram.metapilot.core.opensearch.ExpectedPeakSelector.ExpectedShift;
import uk.ac.quadram.metapilot.task.dda.mag.par.DDAMagOpenSearchParameters;
import uk.ac.quadram.metapilot.task.io.qpx.UnimodCatalog;

/** Stage 1 and the fragment-gain gate take their masses from UNIMOD or the engine, not a hand-picked list. */
class NoHandPickedMassListTest {

	private static UnimodCatalog.Entry byTitle(String title) {
		UnimodCatalog.Entry e = UnimodCatalog.get().byTitle(title);
		assertNotNull(e, title);
		return e;
	}

	@Test
	@DisplayName("UNIMOD classifications are loaded; isotopic labels are recognised")
	void isotopicLabelsAreRecognised() {
		assertTrue(byTitle("Label:13C(4)").isIsotopicLabelOnly());
		assertTrue(byTitle("Label:18O(2)").isIsotopicLabelOnly());
		assertFalse(byTitle("Oxidation").isIsotopicLabelOnly());
		assertFalse(byTitle("Deamidated").isIsotopicLabelOnly());
		long labels = UnimodCatalog.get().entries().stream().filter(UnimodCatalog.Entry::isIsotopicLabelOnly).count();
		assertEquals(126, labels, "records whose every specificity is classed Isotopic label");
	}

	@Test
	@DisplayName("Shift explanations and the isotope-link guard ignore isotopic labels")
	void explanationLookupsIgnoreIsotopicLabels() {
		UnimodCatalog u = UnimodCatalog.get();
		// Label:15N(1), +0.997035, on any residue: no label-free entry within 2 mDa.
		assertTrue(u.allByMass(0.997035, UnimodCatalog.Site.residue('K'), 0.002).isEmpty());
		// Label:18O(2) +4.008491 and Label:13C(4) +4.013419: the nearest label-free mass is further.
		assertTrue(u.distanceToNearestMass(4.0104) > 0.002);
		// A label-free entry is still found.
		assertEquals(0.0, u.distanceToNearestMass(15.994915), 1e-6);
		assertFalse(u.allByMass(0.984016, UnimodCatalog.Site.residue('N'), 0.002).isEmpty());
		ShiftLocaliser s = new ShiftLocaliser(20.0);
		assertTrue(s.explanationsAt("PEPTIDEK", 7, 0.997035).isEmpty());
	}

	@Test
	@DisplayName("Stage-1 candidates are UNIMOD masses without isotopic labels, same masses merged")
	void stageOneCandidatesComeFromUnimod() {
		List<ExpectedShift> shifts = ExpectedPeakSelector.unimodShifts();
		double[] masses = shifts.stream().mapToDouble(s -> s.massDa).toArray();
		for (int i = 1; i < masses.length; i++) {
			assertTrue(masses[i] > masses[i - 1], "sorted and merged");
		}
		ExpectedShift ox = shifts.get(ExpectedPeakSelector.nearest(masses, 15.994915));
		assertEquals(15.994915, ox.massDa, 1e-5);
		assertTrue(ox.name.contains("Oxidation") && ox.name.contains("Phe->Tyr"), ox.name);
		ExpectedShift deam = shifts.get(ExpectedPeakSelector.nearest(masses, 0.984016));
		assertTrue(deam.name.contains("Deamidated"), deam.name);
		// Label:13C(4) (+4.013419) and Label:2H(4) (+4.025107) have no non-label entry at their mass.
		assertTrue(Math.abs(shifts.get(ExpectedPeakSelector.nearest(masses, 4.013419)).massDa - 4.013419) > 1e-3);
		assertTrue(Math.abs(shifts.get(ExpectedPeakSelector.nearest(masses, 4.025107)).massDa - 4.025107) > 1e-3);
	}

	@Test
	void nearestFindsTheClosestSortedValue() {
		double[] m = { -1.0, 0.5, 2.0 };
		assertEquals(0, ExpectedPeakSelector.nearest(m, -5));
		assertEquals(1, ExpectedPeakSelector.nearest(m, 0.6));
		assertEquals(2, ExpectedPeakSelector.nearest(m, 1.4));
		assertEquals(2, ExpectedPeakSelector.nearest(m, 9));
		assertEquals(-1, ExpectedPeakSelector.nearest(new double[0], 1));
	}

	@Test
	@DisplayName("Stage-1 expected shifts do not run after a closed marker search")
	void closedMarkerSearchSkipsExpectedShifts() {
		DDAMagOpenSearchParameters p = new DDAMagOpenSearchParameters();
		p.setRound1SearchMode("closed");
		assertFalse(p.toExpectedPeakConfig().enabled);
		p.setRound1SearchMode("open");
		assertTrue(p.toExpectedPeakConfig().enabled);
	}

	@Test
	@DisplayName("The gate's exemption is the engine's delta_mass_exclude_ranges in the bundled parameters")
	void exemptRangeMatchesBundledParameters() throws IOException {
		for (String resource : new String[] { "/msfragger/open_fragger.params",
				"/msfragger/open_fragger_standalone.params" }) {
			String line = null;
			try (InputStream in = NoHandPickedMassListTest.class.getResourceAsStream(resource)) {
				assertNotNull(in, resource);
				BufferedReader br = new BufferedReader(new InputStreamReader(in, StandardCharsets.UTF_8));
				String l;
				while ((l = br.readLine()) != null) {
					if (l.startsWith("delta_mass_exclude_ranges")) {
						line = l;
					}
				}
			}
			assertNotNull(line, resource);
			String range = line.substring(line.indexOf('(') + 1, line.indexOf(')'));
			String[] lohi = range.split(",");
			assertEquals(DeltaMassEvidenceGate.DELTA_MASS_EXCLUDE_LOWER_DA, Double.parseDouble(lohi[0].trim()), 0.0,
					resource);
			assertEquals(DeltaMassEvidenceGate.DELTA_MASS_EXCLUDE_UPPER_DA, Double.parseDouble(lohi[1].trim()), 0.0,
					resource);
		}
	}
}
