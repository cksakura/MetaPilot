package uk.ac.quadram.metapilot.core.runner;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.io.File;
import java.nio.file.Files;
import java.util.Arrays;
import java.util.List;
import java.util.Map;

import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;

import uk.ac.quadram.metapilot.core.opensearch.OpenSearchLayout;
import uk.ac.quadram.metapilot.core.opensearch.OpenSearchPsm;
import uk.ac.quadram.metapilot.core.opensearch.SubstitutionReassigner;

class SubstitutionQuantPsmTest {

	@TempDir
	File tmp;

	private static final String VALIDATED_HEADER = "File Name\tScan number\tBase Sequence\tFull Sequence\t"
			+ "Peptide Monoisotopic Mass\tScan Retention Time\tPrecursor Charge\tProtein Accession\tScore\tPEP\t"
			+ "Peptide Mass\tMass difference\tMissed cleavages\thit_rank\tis_target";

	@Test
	void variantsOutsideValidatedAndTheirUnmodifiedReferencesAreWritten() throws Exception {
		OpenSearchPsm kept = psm("S1.mzML", 10, "PEPTIDEK", true, "PEPTIDEK[+14.0157]");
		kept.setSubstitutionVerdict(SubstitutionReassigner.Verdict.SUBSTITUTION_NOT_IN_DB.name());
		kept.setResolvedMonoisotopicMass(900.1);
		OpenSearchPsm already = psm("S1.mzML", 11, "OTHERK", true, "OTHERK[+14.0157]");
		already.setSubstitutionVerdict(SubstitutionReassigner.Verdict.SUBSTITUTION_IN_DB.name());
		OpenSearchPsm decoy = psm("S1.mzML", 12, "DECOYK", false, "DECOYK[+14.0157]");
		decoy.setSubstitutionVerdict(SubstitutionReassigner.Verdict.SUBSTITUTION_NOT_IN_DB.name());
		OpenSearchPsm plain = psm("S1.mzML", 13, "PLAINR", true, "PLAINR");
		OpenSearchPsm reference = psm("S2.mzML", 20, "PEPTIDEK", true, "PEPTIDEK");
		OpenSearchPsm oxidised = psm("S2.mzML", 21, "PEPTIDEK", true, "PEPTIDEK[+15.9949]");

		File out = new File(tmp, "substitution_psms.tsv");
		SubstitutionQuantification.Counts counts = SubstitutionQuantification.writePsms(
				Arrays.asList(kept, already, decoy, plain), List.of(already, plain, reference, oxidised),
				Map.of(kept, "V7"), out);

		assertEquals(1, counts.variants);
		assertEquals(1, counts.references);
		List<String> lines = Files.readAllLines(out.toPath());
		assertEquals(SubstitutionQuantification.HEADER, lines.get(0));
		assertEquals(3, lines.size());
		assertTrue(lines.get(1).startsWith("S1.mzML\t10\tPEPTIDEK\tPEPTIDEK[+14.0157]\t900.1\t"));
		assertTrue(lines.get(1).endsWith("\tUNASSIGNED\t3.0\t0.01\t914.0\t14.0\t0\tvariant\tV7"), lines.get(1));
		assertTrue(lines.get(2).startsWith("S2.mzML\t20\tPEPTIDEK\tPEPTIDEK\t"));
		assertTrue(lines.get(2).endsWith("\treference\t"), lines.get(2));
		assertFalse(String.join("\n", lines).contains("GENOME1"));
		assertEquals(1, SubstitutionQuantification.variantRows(out));
		assertTrue(SubstitutionQuantification.hasRoles(out));
	}

	@Test
	void reportFallbackUsesUnexplainedRowsAndAddsReferences() throws Exception {
		File forms = write("forms.tsv", "file\tscan\thit_rank\tsubstitution_verdict\tvariant_id\n"
				+ "S1.mzML\t10\t1\tSUBSTITUTION_NOT_IN_DB\tV3\n"
				+ "S1.mzML\t11\t1\tSUBSTITUTION_IN_DB\t\n"
				+ "S1.mzML\t12\t1\t\t\n");
		File validated = write("validated.tsv", VALIDATED_HEADER + "\n"
				+ "S1.mzML\t11\tOTHERK\tOTHERK[+14.0157]\t1\t1\t2\tGENOME1_1\t1\t0.01\t1\t14\t0\t1\ttrue\n"
				+ "S2.mzML\t30\tPEPTIDEK\tPEPTIDEK\t900\t12\t2\tGENOME1_1\t5\t0.001\t900\t0\t0\t1\ttrue\n"
				+ "S2.mzML\t31\tPEPTIDEK\tPEPTIDEK\t900\t12\t2\trev_G\t1\t0.01\t900\t0\t0\t1\tfalse\n");
		File unexplained = write("unexplained.tsv", VALIDATED_HEADER + "\n"
				+ "S1.mzML\t10\tPEPTIDEK\tPEPTIDEK[+14.0157]\t900.1\t12.5\t2\tGENOME1_1\t3\t0.01\t914\t14\t0\t1\ttrue\n"
				+ "S1.mzML\t11\tOTHERK\tOTHERK[+14.0157]\t1\t1\t2\tGENOME1_1\t1\t0.01\t1\t14\t0\t1\ttrue\n"
				+ "S1.mzML\t14\tDECOYK\tDECOYK[+14.0157]\t1\t1\t2\trev_G\t1\t0.01\t1\t14\t0\t1\tfalse\n");

		File out = new File(tmp, "sv/substitution_psms.tsv");
		SubstitutionQuantification.Counts counts = SubstitutionQuantification.writePsmsFromReports(forms,
				unexplained, validated, out);

		assertEquals(1, counts.variants);
		assertEquals(1, counts.references);
		List<String> lines = Files.readAllLines(out.toPath());
		assertTrue(lines.get(1).startsWith("S1.mzML\t10\tPEPTIDEK\tPEPTIDEK[+14.0157]\t900.1\t"));
		assertTrue(lines.get(1).endsWith("\tUNASSIGNED\t3\t0.01\t914\t14\t0\tvariant\tV3"), lines.get(1));
		assertTrue(lines.get(2).startsWith("S2.mzML\t30\tPEPTIDEK\tPEPTIDEK\t900\t"));
		assertTrue(lines.get(2).endsWith("\treference\t"), lines.get(2));
	}

	@Test
	void annotationAddsRolesRatiosAndSharedMainFeatures() throws Exception {
		File psms = write("psms.tsv", SubstitutionQuantification.HEADER + "\n"
				+ "S1\t10\tPEPTIDEK\tPEPTIDEK[+14.0157]\t914\t12\t2\tUNASSIGNED\t3\t0.01\t914\t14\t0\tvariant\tV3\n"
				+ "S1\t11\tPEPTIDEK\tPEPTIDEK[+14.0157]\t914\t12\t2\tUNASSIGNED\t3\t0.01\t914\t14\t0\tvariant\tV4\n"
				+ "S1\t30\tPEPTIDEK\tPEPTIDEK\t900\t12\t2\tUNASSIGNED\t5\t0.001\t900\t0\t0\treference\t\n");
		File peptides = write("QuantifiedPeptides.tsv",
				"Sequence\tBase Sequence\tProtein Groups\tIntensity_S1\tIntensity_S2\tDetection Type_S1\tDetection Type_S2\n"
						+ "PEPTIDEK[+14.0157]\tPEPTIDEK\tUNASSIGNED\t200\t300\tMSMS\tMBR\n"
						+ "PEPTIDEK\tPEPTIDEK\tUNASSIGNED\t1000\t500\tMSMS\tMSMS\n");
		String peakHeader = "File Name\tBase Sequence\tFull Sequence\tPeak intensity\tPeak RT Start\tPeak RT Apex\t"
				+ "Peak RT End\tPeak MZ\tPeak Charge\tPeak Detection Type\n";
		File side = write("side_peaks.tsv", peakHeader
				+ "S1\tPEPTIDEK\tPEPTIDEK[+14.0157]\t200\t11.8\t12.0\t12.2\t458.0000\t2\tMSMS\n"
				+ "S2\tPEPTIDEK\tPEPTIDEK[+14.0157]\t300\t11.8\t12.0\t12.2\t458.0000\t2\tMBR\n"
				+ "S1\tPEPTIDEK\tPEPTIDEK\t1000\t11.8\t12.0\t12.2\t451.0000\t2\tMSMS\n");
		File main = write("main_peaks.tsv", peakHeader
				+ "S1\tOTHERPEPK\tOTHERPEPK\t900\t11.5\t11.9\t12.4\t458.0030\t2\n"
				+ "S1\tFARK\tFARK\t900\t11.5\t11.9\t12.4\t458.0200\t2\n"
				+ "S1\tLATEK\tLATEK\t900\t14.0\t14.2\t14.4\t458.0000\t2\n"
				+ "S1\tCHARGEK\tCHARGEK\t900\t11.5\t11.9\t12.4\t458.0000\t3\tMSMS\n"
				+ "S2\tMBRPEAKK\tMBRPEAKK\t900\t11.5\t11.9\t12.4\t458.0000\t2\tMSMS\n");

		int sharing = SubstitutionQuantification.annotate(psms, peptides, side, main);

		assertEquals(1, sharing);
		List<String> lines = Files.readAllLines(peptides.toPath());
		assertTrue(lines.get(0).endsWith(
				"\tquant_role\tvariant_ids\treference_sequence\tRatio_S1\tRatio_S2\tshares_main_feature"));
		// S2 was a match-between-runs value: dropped, so no S2 ratio and no S2 shared feature.
		assertTrue(lines.get(1).startsWith("PEPTIDEK[+14.0157]\tPEPTIDEK\tUNASSIGNED\t200\t0\tMSMS\tMBR_removed\t"),
				lines.get(1));
		assertTrue(lines.get(1).endsWith("\tvariant\tV3;V4\tPEPTIDEK\t0.200000\t\tS1:OTHERPEPK"), lines.get(1));
		assertTrue(lines.get(2).endsWith("\treference\t\t\t\t\t"), lines.get(2));
	}

	@Test
	void newSubstitutionFormsGoToSequenceVariantsAndAnOldCopyStays() throws Exception {
		File engine = new File(tmp, "opensearch_MAG_search");
		File combined = OpenSearchLayout.combinedDir(engine);
		assertTrue(combined.mkdirs());
		File legacy = new File(combined, "substitution_forms.tsv");
		Files.writeString(legacy.toPath(), "file\tscan\nS1\t1\n");

		assertEquals(legacy, OpenSearchLayout.substitutionForms(engine, null));
		assertTrue(legacy.isFile());
		assertFalse(OpenSearchLayout.substitutionFormsOutput(engine, null).isFile());

		File output = OpenSearchLayout.substitutionFormsOutput(engine, null);
		assertTrue(output.getParentFile().mkdirs());
		Files.writeString(output.toPath(), "file\tscan\nS2\t2\n");
		assertEquals(output, OpenSearchLayout.substitutionForms(engine, null));
		assertEquals("file\tscan\nS1\t1\n", Files.readString(legacy.toPath()).replace("\r\n", "\n"));
	}

	@Test
	void substitutionPeaksAndCandidatesLiveInSequenceVariants() throws Exception {
		File engine = new File(tmp, "opensearch_MAG_search");
		File combined = OpenSearchLayout.combinedDir(engine);
		assertTrue(combined.mkdirs());
		Files.writeString(new File(combined, "substitution_peaks.tsv").toPath(), "old");

		File peaks = OpenSearchLayout.substitutionPeaks(engine, null);
		assertEquals(new File(engine, "sequence_variants/substitution_peaks.tsv"), peaks);
		assertTrue(peaks.isFile(), "an earlier run's copy is moved");
		assertFalse(new File(combined, "substitution_peaks.tsv").exists());
		assertEquals(new File(engine, "sequence_variants/substitution_candidates.fasta"),
				OpenSearchLayout.substitutionCandidates(engine, null));
	}

	@Test
	void spectraFolderGuardRestoresTheDesignAndLeavesOtherFiles() throws Exception {
		File spectra = new File(tmp, "spectra");
		assertTrue(spectra.mkdirs());
		File mzml = new File(spectra, "S1.mzML");
		Files.writeString(mzml.toPath(), "mzML");
		File design = new File(spectra, "ExperimentalDesign.tsv");
		Files.writeString(design.toPath(), "original");
		long modified = design.lastModified();

		SubstitutionQuantification.SpectraFolderGuard guard = SubstitutionQuantification.SpectraFolderGuard
				.take(spectra);
		Files.writeString(design.toPath(), "rewritten by FlashLFQ");
		Files.writeString(new File(spectra, "other_job.txt").toPath(), "new");
		List<String> appeared = guard.restore();

		assertEquals("original", Files.readString(design.toPath()).replace("\r\n", "\n"));
		assertEquals(modified, design.lastModified());
		assertTrue(mzml.isFile());
		assertTrue(new File(spectra, "other_job.txt").isFile(), "another job's file is not deleted");
		assertEquals(List.of("other_job.txt"), appeared);
	}

	@Test
	void spectraFolderGuardRemovesADesignItCreated() throws Exception {
		File spectra = new File(tmp, "spectra2");
		assertTrue(spectra.mkdirs());
		SubstitutionQuantification.SpectraFolderGuard guard = SubstitutionQuantification.SpectraFolderGuard
				.take(spectra);
		Files.writeString(new File(spectra, "ExperimentalDesign.tsv").toPath(), "written by the side run");
		assertTrue(guard.restore().isEmpty());
		assertFalse(new File(spectra, "ExperimentalDesign.tsv").exists());
	}

	private File write(String name, String content) throws Exception {
		File f = new File(tmp, name);
		Files.writeString(f.toPath(), content);
		return f;
	}

	private static OpenSearchPsm psm(String file, int scan, String peptide, boolean target, String fullSequence) {
		OpenSearchPsm psm = new OpenSearchPsm();
		psm.setFileName(file);
		psm.setScan(scan);
		psm.setHitRank(1);
		psm.setPeptide(peptide);
		psm.setTarget(target);
		psm.setProtein("GENOME1_1");
		psm.setCharge(2);
		psm.setRetentionTime(120.0);
		psm.setPrecursorNeutralMass(914.0);
		psm.setMassDiff(14.0);
		psm.setCalcNeutralPepMass(900.0);
		psm.setQValue(0.01);
		psm.setClassScore(3.0);
		psm.setResolvedModifiedSequence(fullSequence);
		return psm;
	}
}
