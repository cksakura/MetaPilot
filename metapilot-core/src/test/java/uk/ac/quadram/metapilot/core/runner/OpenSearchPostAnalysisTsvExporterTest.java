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
import uk.ac.quadram.metapilot.core.opensearch.ShiftLocaliser;
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
		// Peptide letters stay upper case (OpenSearchPsm#getModificationFamilyKey): the localization
		// guess never enters a peptide's identity; before resolution the bracket carries no position.
		assertEquals("PEPTIDE[+15.9949]", OpenSearchPostAnalysisTsvExporter.canonicalModifiedSequenceOf(mono));
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

	private static final double A_TO_T = 101.04768 - 71.03711;

	private static OpenSearchPsm localised(double massDiff, String msfraggerSites, double classScore, int site) {
		OpenSearchPsm psm = new OpenSearchPsm();
		psm.setPeptide("SAMPLEK");
		psm.setTarget(true);
		psm.setMassDiff(massDiff);
		psm.setNormalizedMassDiff(massDiff);
		psm.setPeakId(-1);
		psm.setCalcNeutralPepMass(700.0);
		psm.setPrecursorNeutralMass(700.0 + massDiff);
		psm.setLocalizedSequence(msfraggerSites);
		psm.setClassScore(classScore);
		psm.setShiftLocalisation(new ShiftLocaliser.Result(8.0, site, 0.97, ShiftLocaliser.State.LOCALISED_EXPLAINED,
				List.of("subst A->T"), A_TO_T));
		return psm;
	}

	@Test
	void localisedPsmsCarryTheLocalisedSiteAndExplainedMass() {
		// MSFragger placed the shift on S (residue 1); the fragments placed it on A (residue 2).
		OpenSearchPsm best = localised(30.0112, "sAMPLEK", 0.9, 1);
		OpenSearchPsm other = localised(30.0098, "SaMPLEK", 0.5, 1);
		assertEquals("SAMPLEK[+30.0106]", best.getModificationFamilyKey(), "one family, whatever the observed mass");
		assertEquals(best.getModificationFamilyKey(), other.getModificationFamilyKey());

		OpenSearchPostAnalysisTsvExporter.resolveCanonicalModifiedSequences(Arrays.asList(best, other));
		for (OpenSearchPsm psm : Arrays.asList(best, other)) {
			assertEquals("SA[+30.0106]MPLEK", psm.getCanonicalModifiedSequence());
			assertEquals(1, psm.getResolvedSiteIndex());
			assertEquals(700.0 + A_TO_T, psm.getFlashLfqMonoisotopicMass(), 1e-9);
			assertEquals("subst A->T", psm.getResolvedShiftExplanation());
		}
		Map<String, String> explanations = OpenSearchPostAnalysisTsvExporter.shiftExplanationsBySequence(
				Arrays.asList(best, other));
		assertEquals(Map.of("SA[+30.0106]MPLEK", "subst A->T"), explanations);
		assertEquals("localised_explained\t2\t0.9700\tsubst A->T\tsubstitution",
				MagOpenSearchHeadlessRunner.shiftColumns(best));
	}

	@Test
	void unlocalisedShiftsKeepTheUnexplainedTag() {
		OpenSearchPsm psm = localised(30.0112, "sAMPLEK", 0.9, 1);
		psm.setShiftLocalisation(new ShiftLocaliser.Result(0.4, 1, 0.5, ShiftLocaliser.State.AMBIGUOUS, List.of(),
				Double.NaN));
		assertEquals("SAMPLEK[unexplained:+30.0112]", psm.getModificationFamilyKey());
		assertEquals("ambiguous\t\t\t\t", MagOpenSearchHeadlessRunner.shiftColumns(psm));
		OpenSearchPostAnalysisTsvExporter.resolveCanonicalModifiedSequences(Arrays.asList(psm));
		assertEquals("", psm.getResolvedShiftExplanation());
		assertEquals(Map.of(), OpenSearchPostAnalysisTsvExporter.shiftExplanationsBySequence(Arrays.asList(psm)));
	}

	private static OpenSearchPsm acceptedShift(ShiftLocaliser.State state, int site) {
		OpenSearchPsm psm = new OpenSearchPsm();
		psm.setPeptide("AAAAGELK");
		psm.setTarget(true);
		psm.setMassDiff(1.9729);
		psm.setNormalizedMassDiff(0.9695);
		psm.setPeakIsotopeOffset(1);
		psm.setPeakId(7);
		psm.setPeakCentroid(0.9841);
		psm.setCalcNeutralPepMass(700.0);
		psm.setPrecursorNeutralMass(701.9729);
		psm.setLocalizedSequence("AAAAGELK");
		psm.setClassScore(0.9);
		psm.setIdentificationGroup("OPEN_ACCEPTED_SHIFT");
		psm.setShiftLocalisation(new ShiftLocaliser.Result(state == ShiftLocaliser.State.UNSUPPORTED ? -3.0 : 6.0,
				site, 0.9, state, List.of(), Double.NaN));
		return psm;
	}

	@Test
	void unsupportedAcceptedShiftsAreExportedUnmodified() {
		OpenSearchPsm psm = acceptedShift(ShiftLocaliser.State.UNSUPPORTED, 7);
		assertEquals("AAAAGELK", psm.getModificationFamilyKey(), "no '+0.9841' claim without fragment support");
		OpenSearchPostAnalysisTsvExporter.resolveCanonicalModifiedSequences(Arrays.asList(psm));
		assertEquals("AAAAGELK", psm.getCanonicalModifiedSequence());
		assertEquals(-1, psm.getResolvedSiteIndex());
		assertEquals(700.0, psm.getFlashLfqMonoisotopicMass(), 1e-9, "FlashLFQ looks for the unmodified peptide");
		assertEquals("unsupported\t\t\t\t", MagOpenSearchHeadlessRunner.shiftColumns(psm));
		assertEquals(1.9729, psm.getMassDiff(), 0.0, "the observed shift stays recorded");
	}

	@Test
	void contradictedAcceptedShiftsSitOnTheirOwnSiteWithoutExplanation() {
		OpenSearchPsm psm = acceptedShift(ShiftLocaliser.State.LOCALISED_UNEXPLAINED, 2);
		OpenSearchPostAnalysisTsvExporter.resolveCanonicalModifiedSequences(Arrays.asList(psm));
		assertEquals("AAA[+0.9841]AGELK", psm.getCanonicalModifiedSequence(), "the fragment-supported site, peak mass");
		assertEquals("", psm.getResolvedShiftExplanation());
		assertEquals("localised_unexplained\t3\t0.9000\t\t", MagOpenSearchHeadlessRunner.shiftColumns(psm));
	}

	@Test
	void peptideReportAddsTheExplanationColumnOnlyWhenAsked(@TempDir Path dir) throws Exception {
		File plain = dir.resolve("plain.tsv").toFile();
		File withShift = dir.resolve("shift.tsv").toFile();
		DdaMagReportExportService.writePeptideReport(plain, null, null, null, null, null, null, "label-free", null,
				null);
		DdaMagReportExportService.writePeptideReport(withShift, null, null, null, null, null, null, "label-free",
				null, null, OpenSearchPostAnalysisTsvExporter.SHIFT_EXPLANATION_COLUMN, Map.of());
		String header = Files.readAllLines(plain.toPath()).get(0);
		assertEquals(header + "\tShift explanation", Files.readAllLines(withShift.toPath()).get(0));
	}

	private static OpenSearchPsm psm(boolean target) {
		OpenSearchPsm psm = new OpenSearchPsm();
		psm.setPeptide("PEPTIDE");
		psm.setTarget(target);
		return psm;
	}
}
