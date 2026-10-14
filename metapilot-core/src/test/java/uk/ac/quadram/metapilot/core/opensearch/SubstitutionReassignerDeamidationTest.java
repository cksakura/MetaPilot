package uk.ac.quadram.metapilot.core.opensearch;

import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.util.ArrayList;
import java.util.List;
import java.util.Set;

import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Test;

/**
 * N-&gt;D and Q-&gt;E are the same molecule as deamidated N and Q, so a
 * deamidation peak is never withdrawn as a substitution.
 */
class SubstitutionReassignerDeamidationTest {

	private static void addPeak(List<OpenSearchPsm> psms, int peakId, double centroid, String peptide,
			String localized) {
		for (int i = 0; i < 60; i++) {
			OpenSearchPsm psm = new OpenSearchPsm();
			psm.setPeakId(peakId);
			psm.setPeakCentroid(centroid);
			psm.setPeptide(peptide);
			psm.setLocalizedSequence(localized);
			psm.setTarget(true);
			psms.add(psm);
		}
	}

	private static SubstitutionReassigner.Config withoutDatabaseTest() {
		SubstitutionReassigner.Config cfg = new SubstitutionReassigner.Config();
		cfg.minInDatabaseFraction = 0.0;
		return cfg;
	}

	@Test
	@DisplayName("Deamidation localized on N is a modification, not N->D")
	void deamidationIsNeverWithdrawn() {
		List<OpenSearchPsm> psms = new ArrayList<>();
		addPeak(psms, 1, 0.98402, "PEPNTIDEK", "PEPnTIDEK");
		addPeak(psms, 2, 14.01565, "PEPGTIDEK", "PEPgTIDEK");

		SubstitutionReassigner.Result result = SubstitutionReassigner.analyse(psms, null, withoutDatabaseTest(),
				null);

		assertFalse(result.substitutionPeakIds.contains(1));
		assertTrue(result.substitutionPeakIds.contains(2), "G->A at +14.0157 is still tested");
	}

	@Test
	@DisplayName("Clearing modificationEquivalent restores the substitution call on the deamidation peak")
	void equivalenceCanBeDisabled() {
		List<OpenSearchPsm> psms = new ArrayList<>();
		addPeak(psms, 1, 0.98402, "PEPNTIDEK", "PEPnTIDEK");
		SubstitutionReassigner.Config cfg = withoutDatabaseTest();
		cfg.modificationEquivalent = Set.of();

		SubstitutionReassigner.Result result = SubstitutionReassigner.analyse(psms, null, cfg, null);

		assertTrue(result.substitutionPeakIds.contains(1));
	}
}
