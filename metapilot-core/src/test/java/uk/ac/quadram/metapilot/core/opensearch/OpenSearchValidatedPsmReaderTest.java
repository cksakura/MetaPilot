package uk.ac.quadram.metapilot.core.opensearch;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNotNull;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.nio.file.Files;
import java.nio.file.Path;
import java.util.List;

import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;

class OpenSearchValidatedPsmReaderTest {

	@Test
	void readsValidatedRowAndKeepsFullSequence(@TempDir Path dir) throws Exception {
		Path tsv = dir.resolve("validated_psms.tsv");
		Files.write(tsv, List.of(
				"File Name\tScan number\tBase Sequence\tFull Sequence\tPeptide Monoisotopic Mass\t"
						+ "Scan Retention Time\tPrecursor Charge\tProtein Accession\tScore\tPEP\t"
						+ "Peptide Mass\tMass difference\tMissed cleavages\thit_rank\tengine\t"
						+ "engine_support\tlocalized_sequence\tdelegate_protein\tis_target\t"
						+ "normalized_mass_diff\tobserved_peak_id\tpeak_id\tpeak_isotope_offset\t"
						+ "hyperscore\texpect\tdelta_mass_score_gain\tlocalization_delta_score\t"
						+ "class_score\tclass_logit\tq_value\tq_value_scope\tvalidation_tier\tstage\t"
						+ "id_group",
				"YS-194-1-1\t12\tPEPTIDEK\tPEP[+15.9949]TIDEK\t900.1\t1.5\t2\t"
						+ "MGYG000000001_00001;MGYG000000002_00001\t0.99\t0.001\t900.2\t15.9949\t0\t1\t"
						+ "MSFRAGGER\tMSFRAGGER\t\tMGYG000000001_00001\ttrue\t15.9949\t1\t1\t0\t"
						+ "20\t1e-5\t\t\t0.99\t2.0\t0.001\tJOINT_GROUP_OPEN_LOCALISED\tPEAK_VALIDATED\t"
						+ "open\tOPEN_LOCALISED"));
		OpenSearchValidatedPsmReader.Result loaded = OpenSearchValidatedPsmReader.read(tsv.toFile());
		assertEquals(1, loaded.psms.size());
		OpenSearchPsm psm = loaded.psms.get(0);
		assertEquals("YS-194-1-1", psm.getFileName());
		assertEquals(12, psm.getScan());
		assertEquals("PEPTIDEK", psm.getPeptide());
		assertEquals("PEP[+15.9949]TIDEK", psm.getCanonicalModifiedSequence());
		assertEquals(900.1, psm.getFlashLfqMonoisotopicMass(), 1e-9);
		assertEquals(90.0, psm.getRetentionTime(), 1e-9);
		assertEquals(2, psm.getAllProteins().length);
		assertEquals("OPEN_LOCALISED", psm.getIdentificationGroup());
		assertTrue(loaded.identificationOnly.isEmpty());
	}

	@Test
	void overlayFillsShiftFromScoredWinners(@TempDir Path dir) throws Exception {
		Path validated = dir.resolve("validated_psms.tsv");
		Files.write(validated, List.of(
				"File Name\tScan number\tBase Sequence\tFull Sequence\tPeptide Monoisotopic Mass\t"
						+ "Scan Retention Time\tPrecursor Charge\tProtein Accession\tScore\tPEP\t"
						+ "Peptide Mass\tMass difference\tMissed cleavages\thit_rank\tengine\t"
						+ "engine_support\tlocalized_sequence\tdelegate_protein\tis_target\t"
						+ "normalized_mass_diff\tobserved_peak_id\tpeak_id\tpeak_isotope_offset\t"
						+ "hyperscore\texpect\tdelta_mass_score_gain\tlocalization_delta_score\t"
						+ "class_score\tclass_logit\tq_value\tq_value_scope\tvalidation_tier\tstage\t"
						+ "id_group",
				"S1\t5\tAAAAK\tA[+40.0313]AAAK\t500.2\t0.1\t2\tMGYG000000001_00001\t0.9\t0.01\t"
						+ "500.3\t40.0313\t0\t1\tMSFRAGGER\tMSFRAGGER\t\tMGYG000000001_00001\ttrue\t"
						+ "40.0313\t-1\t-1\t0\t10\t1e-3\t\t\t0.9\t1.0\t0.01\tJOINT_GROUP_OPEN_LOCALISED\t"
						+ "PEAK_VALIDATED\topen\tOPEN_LOCALISED"));
		Path scored = dir.resolve("scored_psms.tsv");
		Files.write(scored, List.of(
				"file\tscan\tpeptide\tprotein\tis_target\tcharge\tprecursor_neutral_mass\t"
						+ "retention_time_sec\tfound_by\tid_group\tclosed_hit_rank\topen_hit_rank\t"
						+ "sage_hit_rank\tmass_diff\tnormalized_mass_diff\tpeak_id\t"
						+ "peak_hypothesis_accepted\tclosed_hyperscore\tclosed_expect\topen_hyperscore\t"
						+ "open_expect\tdelta_mass_score_gain\tlocalized_sequence\tclass_score\t"
						+ "class_logit\tis_spectrum_winner\tspectrum_q_value\tpeptide_q_value\tvalidated\t"
						+ "stage\tms2_entropy\tms2_matched_fraction\tms2_entropy_margin\trt_error_norm\t"
						+ "query_reassigned\treset_role\tshift_state\tshift_site\tshift_site_probability\t"
						+ "shift_support\tshift_explanation",
				"S1\t5\tAAAAK\tMGYG000000001_00001\ttrue\t2\t500.3\t6\tMSFRAGGER\tOPEN_LOCALISED\t"
						+ "\t1\t\t40.0313\t40.0313\t-1\tfalse\t\t\t10\t1e-3\t\t\t0.9\t1.0\ttrue\t0.01\t"
						+ "0.01\ttrue\topen\t\t\t\t\tfalse\testimating\tlocalised_explained\t1\t0.91\t4.2\t"
						+ "subst G->P|unimod Gly->Pro"));
		OpenSearchValidatedPsmReader.Result loaded = OpenSearchValidatedPsmReader.read(validated.toFile());
		int n = OpenSearchValidatedPsmReader.overlayShifts(loaded.psms, scored.toFile());
		assertEquals(1, n);
		OpenSearchPsm psm = loaded.psms.get(0);
		assertNotNull(psm.getShiftLocalisation());
		assertEquals(ShiftLocaliser.State.LOCALISED_EXPLAINED, psm.getShiftLocalisation().state);
		assertEquals(0, psm.getShiftLocalisation().site);
		assertEquals("subst G->P|unimod Gly->Pro", psm.getResolvedShiftExplanation());
		assertEquals("A[+40.0313]AAAK", psm.getCanonicalModifiedSequence());
	}
}
