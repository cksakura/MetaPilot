package uk.ac.quadram.metapilot.core.mag;

import static org.junit.jupiter.api.Assertions.assertArrayEquals;
import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

import org.dom4j.Element;
import org.junit.jupiter.api.Test;

import uk.ac.quadram.metapilot.task.io.pep.MetaPSM;
import uk.ac.quadram.metapilot.task.io.pep.MetaPeptide;

class DdaMagPostAnalysisServiceTest {

	private static final String GENOME = "MGYG000000001";

	/**
	 * Six peptides unique to one genome, and one target peptide that also maps
	 * to a rev_ protein and outscores them. That peptide is target evidence, not
	 * part of the null; the background decoys alone decide the genome q-value and
	 * never reach the peptide reports.
	 */
	@Test
	void genomeBackgroundDecoysDecideGenomeScoringOnly() {
		List<MetaPeptide> peptides = new ArrayList<>();
		List<MetaPSM> psms = new ArrayList<>();
		for (int i = 0; i < 6; i++) {
			add(peptides, psms, "TARGET" + (char) ('A' + i) + "K", 10.0, GENOME + "_0000" + i);
		}
		add(peptides, psms, "SHAREDK", 12.0, GENOME + "_00009", "rev_MGYG000000002_00001");

		DdaMagPostAnalysisService.Result weakNull = analyzeWithDecoys(peptides, psms, 1.0);
		DdaMagPostAnalysisService.Result strongNull = analyzeWithDecoys(peptides, psms, 15.0);

		assertTrue(q(weakNull) <= DdaMagProteinInferenceService.DEFAULT_GENOME_Q_VALUE_THRESHOLD,
				"low-scoring decoys accept the genome");
		assertTrue(q(strongNull) > DdaMagProteinInferenceService.DEFAULT_GENOME_Q_VALUE_THRESHOLD,
				"decoys that outscore every target reject it");

		assertEquals(weakNull.peptideList.size(), strongNull.peptideList.size());
		assertEquals(weakNull.psmMap.keySet(), strongNull.psmMap.keySet());
		assertEquals(weakNull.pepPsmCountMap.keySet(), strongNull.pepPsmCountMap.keySet());
		assertArrayEquals(weakNull.filePsmCount, strongNull.filePsmCount);
		for (MetaPeptide peptide : strongNull.peptideList) {
			assertTrue(!peptide.getSequence().startsWith("DECOY"), peptide.getSequence());
		}
	}

	private static DdaMagPostAnalysisService.Result analyzeWithDecoys(List<MetaPeptide> peptides,
			List<MetaPSM> psms, double decoyScore) {
		List<MetaPeptide> decoyPeptides = new ArrayList<>();
		List<MetaPSM> decoyPsms = new ArrayList<>();
		for (int i = 0; i < 200; i++) {
			add(decoyPeptides, decoyPsms, "DECOY" + i + "K", decoyScore, "rev_MGYG000000003_" + i);
		}
		return DdaMagPostAnalysisService.analyze(input(peptides, psms, decoyPeptides, decoyPsms));
	}

	private static double q(DdaMagPostAnalysisService.Result result) {
		return result.genomeDiaStats.get(GENOME)[DdaMagGenomeScoringService.IDX_Q_VALUE];
	}

	private static DdaMagPostAnalysisService.Input input(List<MetaPeptide> peptides, List<MetaPSM> psms,
			List<MetaPeptide> decoyPeptides, List<MetaPSM> decoyPsms) {
		HashMap<String, Integer> fileIdMap = new HashMap<>();
		fileIdMap.put("S1", 0);
		ArrayList<String> fileList = new ArrayList<>(List.of("S1"));
		return new DdaMagPostAnalysisService.Input(peptides.toArray(new MetaPeptide[0]),
				psms.toArray(new MetaPSM[0]), fileIdMap, fileList, "rev_", "MGYG", null, false,
				decoyPeptides == null ? null : decoyPeptides.toArray(new MetaPeptide[0]),
				decoyPsms == null ? null : decoyPsms.toArray(new MetaPSM[0]));
	}

	private static void add(List<MetaPeptide> peptides, List<MetaPSM> psms, String sequence, double score,
			String... proteins) {
		peptides.add(new MetaPeptide(sequence, sequence, proteins, new double[] { 1000.0 }, null) {
			@Override
			public Object[] getTableObjects() {
				return null;
			}

			@Override
			public Element getXmlPepElement() {
				return null;
			}
		});
		psms.add(new MetaPSM("S1", psms.size() + 1, sequence, sequence, 2, proteins, score, 0.001) {
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
		});
	}
}
