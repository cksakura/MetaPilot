package uk.ac.quadram.metapilot.core.mag;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.util.HashMap;

import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Test;

import uk.ac.quadram.metapilot.task.io.pep.MetaPSM;
import uk.ac.quadram.metapilot.task.io.pep.MetaPeptide;

/**
 * Decoy definition used by the genome and protein resampling nulls in
 * {@link DdaMagGenomeScoringService}.
 */
class DdaMagGenomeScoringDecoyTest {

	private static MetaPSM psm(int scan, String seq, double score, String... proteins) {
		return new MetaPSM("run", scan, seq, seq, 2, proteins, score, 0.0) {
			@Override
			public double getExpMass() {
				return 0;
			}

			@Override
			public double getRt() {
				return 0;
			}

			@Override
			public double getCalMass() {
				return 0;
			}

			@Override
			public double getMassDiff() {
				return 0;
			}

			@Override
			public int getMiss() {
				return 0;
			}
		};
	}

	@Test
	@DisplayName("A peptide is decoy evidence only when all of its proteins are decoys")
	void decoyOnlyNeedsEveryProteinToBeDecoy() {
		assertTrue(DdaMagGenomeScoringService.isDecoyOnly(new String[] { "rev_MGYG1_1", "rev_MGYG2_7" }, "rev_"));
		assertFalse(DdaMagGenomeScoringService.isDecoyOnly(new String[] { "MGYG1_1", "rev_MGYG2_7" }, "rev_"));
		assertFalse(DdaMagGenomeScoringService.isDecoyOnly(new String[] { "MGYG1_1" }, "rev_"));
		assertFalse(DdaMagGenomeScoringService.isDecoyOnly(new String[0], "rev_"));
		assertFalse(DdaMagGenomeScoringService.isDecoyOnly(null, "rev_"));
	}

	@Test
	@DisplayName("A high-scoring target that also maps to a reversed protein stays out of the genome null")
	void sharedTargetDecoySequenceDoesNotEnterGenomeNull() {
		MetaPSM[] psms = {
				psm(1, "PEPTIDEK", 20.0, "MGYG000000001_00001"),
				psm(2, "ARIPLAK", 30.0, "MGYG000000002_00001", "rev_MGYG000000003_00001"),
				psm(3, "DECOYPEPR", 1.0, "rev_MGYG000000004_00001") };

		DdaMagGenomeScoringService.Result result = DdaMagGenomeScoringService.refineGenomes(new MetaPeptide[0],
				"rev_", psms, "MGYG");

		assertEquals(1, result.decoyPeptideCount);
		double[] stats = result.diaStyleGenomeStats.get("MGYG000000001");
		assertEquals(1.0 / 10001.0, stats[DdaMagGenomeScoringService.IDX_P_EMP], 1e-12);
	}

	@Test
	@DisplayName("Protein null uses the same decoy definition")
	void sharedTargetDecoySequenceDoesNotEnterProteinNull() {
		MetaPSM[] psms = {
				psm(1, "PEPTIDEK", 20.0, "MGYG000000001_00001"),
				psm(2, "ARIPLAK", 30.0, "MGYG000000002_00001", "rev_MGYG000000003_00001"),
				psm(3, "DECOYPEPR", 1.0, "rev_MGYG000000004_00001") };

		HashMap<String, double[]> proteins = DdaMagGenomeScoringService.computeProteinScoreMap(new MetaPeptide[0],
				psms, "rev_");

		assertEquals(1.0 / 10001.0, proteins.get("MGYG000000001_00001")[DdaMagGenomeScoringService.IDX_P_EMP],
				1e-12);
		assertTrue(proteins.containsKey("MGYG000000002_00001"));
	}
}
