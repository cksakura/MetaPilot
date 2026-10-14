package uk.ac.quadram.metapilot.core.opensearch;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.io.File;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Map;

import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;

import uk.ac.quadram.metapilot.task.io.pep.MetaPSM;
import uk.ac.quadram.metapilot.task.io.qpx.QpxPsmParquetWriter;

class OpenSearchQpxExportTest {

	@TempDir
	Path tmp;

	@Test
	void writesOpenSearchPsms() throws Exception {
		Map<String, Double> staticMods = Collections.singletonMap("C", 57.02146);
		List<MetaPSM> psms = new ArrayList<>();

		OpenSearchPsm plain = psm("PEPTCIDEK", 1);
		plain.setPeakId(0);
		plain.setPeakCentroid(0.0);
		plain.setNormalizedMassDiff(0.0);
		psms.add(new OpenSearchMetaPsm(plain, staticMods));

		OpenSearchPsm localized = psm("PEPTMIDEK", 2);
		localized.setMassDiff(15.9949);
		localized.setNormalizedMassDiff(15.9949);
		localized.setPeakId(3);
		localized.setPeakCentroid(15.9949);
		localized.setLocalizedSequence("PEPTmIDEK");
		localized.setResolvedSiteIndex(4);
		psms.add(new OpenSearchMetaPsm(localized, staticMods));

		OpenSearchPsm unexplained = psm("PEPTIDEK", 3);
		unexplained.setMassDiff(123.4567);
		unexplained.setPeakId(-1);
		psms.add(new OpenSearchMetaPsm(unexplained, staticMods));

		File out = new File(tmp.toFile(), QpxPsmParquetWriter.DEFAULT_FILE_NAME);
		QpxPsmParquetWriter.Summary summary = QpxPsmParquetWriter.write(out, psms);

		assertEquals(3, summary.written, summary.toString());
		assertTrue(out.length() > 0);
	}

	private static OpenSearchPsm psm(String peptide, int scan) {
		OpenSearchPsm psm = new OpenSearchPsm();
		psm.setPeptide(peptide);
		psm.setTarget(true);
		psm.setFileName("sample_A");
		psm.setScan(scan);
		psm.setCharge(2);
		psm.setCalcNeutralPepMass(1000.0 + scan);
		psm.setPrecursorNeutralMass(1000.0 + scan);
		psm.setRetentionTime(600.0);
		psm.setHyperscore(30.0);
		psm.setExpect(1e-5);
		psm.setClassScore(0.99);
		psm.setClassLogit(4.6);
		psm.setQValue(0.001);
		psm.setProtein("MGYG000000001_00001");
		psm.setEngine(OpenSearchPsm.SearchEngine.MSFRAGGER);
		return psm;
	}
}
