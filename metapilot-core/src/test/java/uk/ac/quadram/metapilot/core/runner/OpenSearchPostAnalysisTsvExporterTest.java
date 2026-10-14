package uk.ac.quadram.metapilot.core.runner;

import static org.junit.jupiter.api.Assertions.assertEquals;

import java.io.File;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Arrays;
import java.util.List;
import java.util.Map;

import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;

import uk.ac.quadram.metapilot.core.aminoacid.AminoAcidProperty;
import uk.ac.quadram.metapilot.core.opensearch.OpenSearchPsm;
import uk.ac.quadram.metapilot.quant.flashLFQ.FlashLfqQuanPeptide;

class OpenSearchPostAnalysisTsvExporterTest {

	@Test
	void postAnalysisReceivesTargetsOnly() {
		OpenSearchPsm target = psm(true);
		OpenSearchPsm decoy = psm(false);

		List<OpenSearchPsm> result = OpenSearchPostAnalysisTsvExporter.targetsOnly(Arrays.asList(target, decoy));

		assertEquals(1, result.size());
		assertEquals(target, result.get(0));
	}

	@Test
	void genomeScoringBackgroundReceivesDecoysOnly() {
		OpenSearchPsm target = psm(true);
		OpenSearchPsm decoy = psm(false);

		List<OpenSearchPsm> result = OpenSearchPostAnalysisTsvExporter.decoysOnly(Arrays.asList(target, decoy));

		assertEquals(1, result.size());
		assertEquals(decoy, result.get(0));
	}

	@Test
	void flashLfqDecoyRowsKeepTheirIntensityForGenomeScoring(@TempDir Path dir) throws Exception {
		File tsv = dir.resolve("QuantifiedPeptides.tsv").toFile();
		Files.write(tsv.toPath(), List.of(
				"Sequence\tBase Sequence\tProtein Groups\tIntensity_S1\tDetection Type_S1",
				"PEPTIDEK\tPEPTIDEK\tMGYG000000001_00001\t1000\tMSMS",
				"KEDITPEP\tKEDITPEP\trev_MGYG000000001_00001\t500\tMSMS"));

		Map<String, FlashLfqQuanPeptide> quantified = OpenSearchPostAnalysisTsvExporter.readQuantifiedPeptides(tsv,
				List.of("S1"));

		assertEquals(1000.0, quantified.get("PEPTIDEK").getIntensity()[0]);
		assertEquals(500.0, quantified.get("KEDITPEP").getIntensity()[0]);
	}

	@Test
	void canonicalKeyUsesFamilyRootMassRatherThanRawPsmMass() {
		OpenSearchPsm mono = psm(true);
		mono.setLocalizedSequence("PEpTIDE");
		mono.setMassDiff(15.9961);
		mono.setPeakId(10);
		mono.setPeakCentroid(15.9949);

		OpenSearchPsm isotope = psm(true);
		isotope.setLocalizedSequence("PEpTIDE");
		isotope.setMassDiff(16.9998);
		isotope.setObservedPeakId(11);
		isotope.setPeakId(10);
		isotope.setPeakCentroid(15.9949);
		isotope.setNormalizedMassDiff(15.9964);
		isotope.setPeakIsotopeOffset(1);

		assertEquals(OpenSearchPostAnalysisTsvExporter.canonicalModifiedSequenceOf(mono),
				OpenSearchPostAnalysisTsvExporter.canonicalModifiedSequenceOf(isotope));
		assertEquals("PEpTIDE[+15.9949]", OpenSearchPostAnalysisTsvExporter.canonicalModifiedSequenceOf(mono));
	}

	@Test
	void isotopeMemberOfZeroFamilyUsesUnmodifiedKey() {
		OpenSearchPsm isotope = psm(true);
		isotope.setLocalizedSequence("pEPTIDE");
		isotope.setMassDiff(AminoAcidProperty.isotope_deltaM);
		isotope.setObservedPeakId(3);
		isotope.setPeakId(2);
		isotope.setPeakCentroid(0.0);
		isotope.setNormalizedMassDiff(0.0);
		isotope.setPeakIsotopeOffset(1);

		assertEquals("PEPTIDE", OpenSearchPostAnalysisTsvExporter.canonicalModifiedSequenceOf(isotope));
	}

	private static OpenSearchPsm psm(boolean target) {
		OpenSearchPsm psm = new OpenSearchPsm();
		psm.setPeptide("PEPTIDE");
		psm.setTarget(target);
		return psm;
	}
}
