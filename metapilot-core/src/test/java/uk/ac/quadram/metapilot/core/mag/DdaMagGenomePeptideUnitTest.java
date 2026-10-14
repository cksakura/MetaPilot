package uk.ac.quadram.metapilot.core.mag;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.util.ArrayList;
import java.util.List;

import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Test;

import uk.ac.quadram.metapilot.task.io.pep.MetaPSM;
import uk.ac.quadram.metapilot.task.io.pep.MetaPeptide;

/**
 * Peptide unit and decoy-genome calibration of {@link DdaMagGenomeScoringService}.
 */
class DdaMagGenomePeptideUnitTest {

	private static final String GENOME = "MGYG000000001";
	private static final String DECOY_GENOME = "MGYG000000009";

	private static MetaPSM psm(int scan, String base, String full, int charge, double score, String protein) {
		return new MetaPSM("run01.mzML", scan, base, full, charge, new String[] { protein }, score, 0.0) {
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

	private static DdaMagGenomeScoringService.Result score(List<MetaPSM> psms) {
		return DdaMagGenomeScoringService.refineGenomes(new MetaPeptide[0], "rev_", psms.toArray(new MetaPSM[0]),
				"MGYG");
	}

	@Test
	@DisplayName("Charge states and modified forms of one peptide are one scoring unit")
	void peptideFormsCountOnce() {
		List<MetaPSM> psms = new ArrayList<>();
		psms.add(psm(1, "PEPTLDEK", "PEPTLDEK", 2, 20.0, GENOME + "_00001"));
		psms.add(psm(2, "PEPTLDEK", "PEPTLDEK", 3, 20.0, GENOME + "_00001"));
		psms.add(psm(3, "PEPTLDEK", "PEPT[+15.9949]LDEK", 2, 20.0, GENOME + "_00001"));
		psms.add(psm(4, "PEPTIDEK", "PEPTIDEK", 2, 20.0, GENOME + "_00001"));
		psms.add(psm(5, "DECOYPEPR", "DECOYPEPR", 2, 1.0, "rev_" + DECOY_GENOME + "_00001"));

		DdaMagGenomeScoringService.Result result = score(psms);

		assertEquals(1, (int) result.diaStyleRuninfo.genomeTotalPeptidesMap.get(GENOME));
	}

	@Test
	@DisplayName("Reversed proteomes go through the genome statistic for calibration, not selection")
	void decoyGenomesAreCalibratedNotSelected() {
		String letters = "ACDEFGHKMNPQRSTVWY";
		List<MetaPSM> psms = new ArrayList<>();
		for (int i = 0; i < letters.length(); i++) {
			String target = "TARGETPEP" + letters.charAt(i) + "K";
			String decoy = "DECOYPEP" + letters.charAt(i) + "R";
			psms.add(psm(i, target, target, 2, 20.0, GENOME + "_00001"));
			psms.add(psm(100 + i, decoy, decoy, 2, 1.0, "rev_" + DECOY_GENOME + "_00001"));
		}

		DdaMagGenomeScoringService.Result result = score(psms);

		assertEquals(letters.length(), (int) result.diaStyleRuninfo.decoyGenomeScoreNMap.get(DECOY_GENOME));
		assertTrue(result.diaStyleRuninfo.decoyGenomeScorePMap.get(DECOY_GENOME) > 0.01);
		assertTrue(result.diaStyleGenomeStats.get(GENOME)[DdaMagGenomeScoringService.IDX_Q_VALUE] <= 0.01);
		assertEquals(1, result.diaStyleGenomeStats.size());
	}
}
