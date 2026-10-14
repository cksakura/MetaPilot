package uk.ac.quadram.metapilot.core.opensearch;

import static org.junit.jupiter.api.Assertions.assertArrayEquals;
import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertNotNull;
import static org.junit.jupiter.api.Assertions.assertNull;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.util.ArrayList;
import java.util.EnumMap;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.dom4j.DocumentHelper;
import org.dom4j.Element;
import org.junit.jupiter.api.Test;

import uk.ac.quadram.metapilot.core.opensearch.ClosedOpenCascade.IdClass;
import uk.ac.quadram.metapilot.core.opensearch.ModificationValidator.MassClass;
import uk.ac.quadram.metapilot.core.opensearch.ModificationValidator.PeakVerdict;

class ClosedOpenCascadeTest {

	private static int nextScan = 1;

	private static OpenSearchPsm psm(String peptide, boolean target, double score, int peakId) {
		OpenSearchPsm psm = new OpenSearchPsm();
		psm.setFileName("run1");
		psm.setScan(nextScan++);
		psm.setPeptide(peptide);
		psm.setTarget(target);
		psm.setClassScore(score);
		psm.setPeakId(peakId);
		return psm;
	}

	private static PeakVerdict verdict(int peakId, double centroid, MassClass massClass, boolean accepted) {
		return new PeakVerdict(peakId, centroid, massClass, "", 0, 0, 0.0, 0.0, Double.NaN, true, accepted, "test",
				0);
	}

	@Test
	void peptideFilterTightensEveryClassWhenTheUnionExceedsTheFdr() {
		Map<IdClass, List<OpenSearchPsm>> pools = new EnumMap<>(IdClass.class);
		for (IdClass c : IdClass.values()) {
			List<OpenSearchPsm> pool = new ArrayList<>();
			for (int i = 0; i < 100; i++) {
				pool.add(psm("PEPTK" + i, true, 200 - i, 0));
			}
			pool.add(psm("DECOY" + c.name(), false, 50, 0));
			pools.put(c, pool);
		}

		ClosedOpenCascade.PeptideFilter filter = ClosedOpenCascade.filterPeptides(pools, pools, 0.01);

		for (IdClass c : IdClass.values()) {
			assertArrayEquals(new int[] { 100, 1 }, filter.peptidesBefore.get(c), c.label);
			assertArrayEquals(new int[] { 100, 0 }, filter.peptidesAfter.get(c), c.label);
			assertArrayEquals(new int[] { 100, 0 }, filter.psmsKept.get(c), c.label);
		}
		assertEquals(100, filter.unionTargetsBefore);
		assertEquals(3, filter.unionDecoysBefore);
		assertEquals(0, filter.unionDecoysAfter);
		assertEquals(0.9, filter.factor, 1e-12);
		assertEquals(300, filter.kept.size());
		assertTrue(filter.kept.stream().allMatch(OpenSearchPsm::isTarget));
	}

	@Test
	void peptideFilterLeavesTheThresholdAloneWhenTheUnionPasses() {
		List<OpenSearchPsm> pool = new ArrayList<>();
		for (int i = 0; i < 100; i++) {
			pool.add(psm("PEPTK" + i, true, 200 - i, 0));
		}
		pool.add(psm("DECOYK", false, 50, 0));
		Map<IdClass, List<OpenSearchPsm>> pools = new EnumMap<>(IdClass.class);
		pools.put(IdClass.CLOSED, pool);

		ClosedOpenCascade.PeptideFilter filter = ClosedOpenCascade.filterPeptides(pools, pools, 0.01);

		assertEquals(1.0, filter.factor, 0.0);
		assertEquals(101, filter.kept.size());
		assertArrayEquals(new int[] { 0, 0 }, filter.psmsKept.get(IdClass.OPEN_KNOWN_SHIFT));
	}

	@Test
	void peptideUnitsMergeIsoleucineAndLeucine() {
		List<OpenSearchPsm> pool = List.of(psm("PEPTIK", true, 10, 0), psm("PEPTLK", true, 5, 0));
		Map<String, Double> q = ClosedOpenCascade.peptideQValues(pool);
		assertEquals(1, q.size());
	}

	@Test
	void openStageIdentifiesOnlyAcceptedZeroAndKnownPeaks() {
		List<OpenSearchPsm> remaining = new ArrayList<>();
		for (int i = 0; i < 20; i++) {
			remaining.add(psm("ZERO" + i, true, 100 - i, 0));
		}
		for (int i = 0; i < 10; i++) {
			remaining.add(psm("KNOWN" + i, true, 100 - i, 1));
		}
		for (int i = 0; i < 10; i++) {
			remaining.add(psm("SHIFT" + i, true, 100 - i, 2));
		}
		for (int i = 0; i < 5; i++) {
			remaining.add(psm("NOPEAK" + i, true, 100 - i, -1));
		}
		ModificationValidator.Gate gate = new ModificationValidator.Gate(Set.of(0, 1),
				List.of(verdict(0, 0.0, MassClass.ZERO, true), verdict(1, 15.9949, MassClass.KNOWN, true),
						verdict(2, 7.3, MassClass.UNEXPLAINED, false)),
				Set.of());

		ModificationValidator.Result result = ModificationValidator.validateRemaining(remaining, null, gate);

		assertEquals(30, result.validated.size());
		assertTrue(result.validated.stream().allMatch(p -> p.getPeakId() == 0 || p.getPeakId() == 1));
		assertEquals(15, result.unexplainedIdentifications.size());
		assertTrue(result.unexplainedIdentifications.stream().noneMatch(result.validated::contains));
		assertEquals(45, result.diagnostic.size());
		PeakVerdict shifted = result.verdicts.stream().filter(v -> v.peakId == 2).findFirst().orElseThrow();
		assertEquals(0, shifted.psmsRetained);
		assertFalse(shifted.accepted);
	}

	@Test
	void gateRejectsWithoutFoldingWhenAsked() {
		List<OpenSearchPsm> psms = new ArrayList<>();
		for (int i = 0; i < 10; i++) {
			OpenSearchPsm p = psm("ZERO" + i, true, 10, 0);
			p.setPeakCentroid(0.0);
			psms.add(p);
		}
		List<OpenSearchPsm> shifted = new ArrayList<>();
		for (int i = 0; i < 60; i++) {
			OpenSearchPsm p = psm("SHIFT" + i, true, 10, 1);
			p.setPeakCentroid(2.718);
			p.setScoreWithDeltaMass(5.0);
			p.setScoreWithoutDeltaMass(5.0);
			shifted.add(p);
		}
		psms.addAll(shifted);
		DeltaMassEvidenceGate.Config cfg = new DeltaMassEvidenceGate.Config();
		cfg.foldRejected = false;

		DeltaMassEvidenceGate.Result result = DeltaMassEvidenceGate.apply(psms, null, cfg, null);

		assertTrue(result.rejectedPeakIds.contains(1));
		assertEquals(0, result.psmsReassigned);
		for (OpenSearchPsm p : shifted) {
			assertEquals(1, p.getPeakId());
			assertEquals(2.718, p.getPeakCentroid(), 0.0);
		}
	}

	@Test
	void closedPsmsCarryTheirVariableModificationsAndIsotopeOffset() {
		OpenSearchPsm p = psm("PEPMK", true, 10, -1);
		p.setClosedStage(true);
		p.setClosedModifications(new double[] { 42.0106, 0, 0, 0, 15.9949, 0 });
		p.setMassDiff(1.0040);
		p.setCalcNeutralPepMass(700.0);

		ClosedOpenCascade.prepareClosed(List.of(p));

		assertEquals(1, p.getPeakIsotopeOffset());
		assertEquals(1.0040 - ClosedOpenCascade.ISOTOPE_SPACING, p.getNormalizedMassDiff(), 1e-9);
		assertEquals(2.0, p.getExternalFeatures().get(ClosedOpenCascade.VARIABLE_MODIFICATIONS_FEATURE), 0.0);
		assertEquals(2, p.getClosedModificationCount());
		assertEquals(42.0106 + 15.9949, p.getClosedModificationMass(), 1e-9);
		assertEquals("[+42.0106]-PEPM[+15.9949]K", p.getClosedModifiedSequence());
		assertEquals(p.getClosedModifiedSequence(), p.getModificationFamilyKey());
		assertEquals(-1, p.getPeakId());
	}

	@Test
	void readerTakesVariableModificationsAndIgnoresFixedCysteine() {
		Element hit = DocumentHelper.createElement("search_hit");
		Element info = hit.addElement("modification_info");
		info.addAttribute("mod_nterm_mass", "43.0184");
		info.addElement("mod_aminoacid_mass").addAttribute("position", "2").addAttribute("mass", "160.0307");
		info.addElement("mod_aminoacid_mass").addAttribute("position", "4").addAttribute("mass", "147.0354")
				.addAttribute("variable", "15.9949");
		info.addElement("mod_aminoacid_mass").addAttribute("position", "5").addAttribute("mass", "147.0354");

		double[] mods = OpenSearchPepxmlReader.readVariableModifications(hit, "ACDMMK");

		assertNotNull(mods);
		assertEquals(7, mods.length);
		assertEquals(42.0106, mods[0], 1e-3);
		assertEquals(0.0, mods[2], 0.0);
		assertEquals(15.9949, mods[4], 1e-4);
		assertEquals(15.9949, mods[5], 1e-3);
	}

	@Test
	void readerReturnsNullForAnUnmodifiedPeptide() {
		Element hit = DocumentHelper.createElement("search_hit");
		hit.addElement("modification_info").addElement("mod_aminoacid_mass").addAttribute("position", "1")
				.addAttribute("mass", "160.0307");
		assertNull(OpenSearchPepxmlReader.readVariableModifications(hit, "CPEPK"));
	}
}
