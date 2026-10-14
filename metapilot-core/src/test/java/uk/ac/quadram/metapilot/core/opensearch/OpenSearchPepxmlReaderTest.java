package uk.ac.quadram.metapilot.core.opensearch;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.nio.file.Files;
import java.nio.file.Path;
import java.util.List;

import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;

class OpenSearchPepxmlReaderTest {

	@TempDir
	Path tempDir;

	@Test
	void readsCurrentPtmResultLocalizationAttributes() throws Exception {
		Path pepxml = writePepXml("<ptm_result localization_peptide=\"PEpTIDE\" "
				+ "score_best_position=\"25.5\" score_second_best_position=\"20.0\" "
				+ "score_all_unshifted=\"10.5\" ions_best_position=\"9\" ions_all_unshifted=\"4\"/>");

		List<OpenSearchPsm> psms = new OpenSearchPepxmlReader("rev_").read(pepxml.toFile());

		assertEquals(1, psms.size());
		OpenSearchPsm psm = psms.get(0);
		assertEquals("PEpTIDE", psm.getLocalizedSequence());
		assertEquals(25.5, psm.getScoreWithDeltaMass(), 1e-12);
		assertEquals(20.0, psm.getSecondBestScoreWithDeltaMass(), 1e-12);
		assertEquals(10.5, psm.getScoreWithoutDeltaMass(), 1e-12);
		assertEquals(15.0, psm.getDeltaMassScoreGain(), 1e-12);
		assertEquals(5.5, psm.getLocalizationDeltaScore(), 1e-12);
		assertEquals(9.0, psm.getIonsWithDeltaMass(), 1e-12);
		assertEquals(4.0, psm.getIonsWithoutDeltaMass(), 1e-12);
		assertEquals(1, psm.getLocalizationSiteCount());
		assertTrue(OpenSearchPepxmlReader.hasLocalizationEvidence(psms));
	}

	@Test
	void readsLegacyPtmResultLocalizationAttributes() throws Exception {
		Path pepxml = writePepXml("<ptm_result localization_peptide=\"pEPTIDE\" "
				+ "best_score_with_ptm=\"19.4\" second_best_score_with_ptm=\"17.2\" "
				+ "score_without_ptm=\"13.8\" ptm_mass=\"15.9949\"/>");

		OpenSearchPsm psm = new OpenSearchPepxmlReader("rev_").read(pepxml.toFile()).get(0);

		assertEquals("pEPTIDE", psm.getLocalizedSequence());
		assertEquals(19.4, psm.getScoreWithDeltaMass(), 1e-12);
		assertEquals(17.2, psm.getSecondBestScoreWithDeltaMass(), 1e-12);
		assertEquals(13.8, psm.getScoreWithoutDeltaMass(), 1e-12);
		assertTrue(psm.getPtmResultAttributes().containsKey("ptm_mass"));
	}

	private Path writePepXml(String ptmResult) throws Exception {
		String xml = "<?xml version=\"1.0\"?>\n"
				+ "<msms_pipeline_analysis><msms_run_summary base_name=\"sample\">"
				+ "<spectrum_query start_scan=\"42\" assumed_charge=\"2\" precursor_neutral_mass=\"1000.0\" "
				+ "retention_time_sec=\"120.0\"><search_result>"
				+ "<search_hit hit_rank=\"1\" peptide=\"PEPTIDE\" peptide_prev_aa=\"K\" peptide_next_aa=\"R\" "
				+ "calc_neutral_pep_mass=\"984.0\" massdiff=\"15.9949\" num_tol_term=\"2\" "
				+ "num_missed_cleavages=\"0\" num_matched_ions=\"8\" tot_num_ions=\"12\" protein=\"target\">"
				+ "<search_score name=\"hyperscore\" value=\"30.0\"/>" + ptmResult
				+ "</search_hit></search_result></spectrum_query></msms_run_summary></msms_pipeline_analysis>";
		Path file = tempDir.resolve("sample.pepXML");
		Files.writeString(file, xml);
		return file;
	}
}
