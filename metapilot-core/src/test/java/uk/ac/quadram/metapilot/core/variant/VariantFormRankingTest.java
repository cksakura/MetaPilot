package uk.ac.quadram.metapilot.core.variant;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;

import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;

class VariantFormRankingTest {

	@TempDir
	File tmp;

	private static VariantCandidate candidate(VariantCandidate.Engine engine) {
		VariantCandidate c = new VariantCandidate();
		c.engine = engine;
		c.file = "s1";
		c.scan = 3;
		c.charge = 2;
		c.basePeptide = "GASDK";
		c.position = 1;
		c.fromAa = 'G';
		c.toAa = 'A';
		c.variantPeptide = "AASDK";
		c.deltaMass = 71.03711 - 57.02146;
		return c;
	}

	@Test
	@DisplayName("Spectra resolve to <base>.mzML, given directly or beside a raw file")
	void mzmlBySample() throws IOException {
		File direct = new File(tmp, "s1.mzML");
		File raw = new File(tmp, "s2.raw");
		File converted = new File(tmp, "s2.mzML");
		File rawOnly = new File(tmp, "s3.raw");
		for (File f : new File[] { direct, raw, converted, rawOnly }) {
			Files.write(f.toPath(), new byte[0]);
		}
		Map<String, File> m = VariantFormRanking.mzmlBySample(List.of(direct, raw, rawOnly));
		assertEquals(direct, m.get("s1"));
		assertEquals("s2.mzml", m.get("s2").getName().toLowerCase());
		assertFalse(m.containsKey("s3"));
	}

	@Test
	@DisplayName("Off, unconfigured or without DIA-NN: open search NOT_RANKED, pFind untouched, nothing thrown")
	void fallbacks() {
		VariantCandidate open = candidate(VariantCandidate.Engine.OPEN_SEARCH);
		VariantCandidate pfind = candidate(VariantCandidate.Engine.PFIND);
		List<VariantCandidate> all = List.of(open, pfind);
		List<String> warnings = new ArrayList<>();

		VariantFormRanking.run(all, null, tmp, s -> {
		}, warnings::add);
		assertEquals(VariantCandidate.Ms2Verdict.NOT_RANKED, open.ms2Verdict);
		assertEquals(VariantCandidate.Ms2Verdict.NOT_APPLICABLE, pfind.ms2Verdict);

		VariantFormRanking.Inputs off = new VariantFormRanking.Inputs();
		off.enabled = false;
		VariantFormRanking.run(all, off, tmp, s -> {
		}, warnings::add);
		assertEquals(VariantCandidate.Ms2Verdict.NOT_RANKED, open.ms2Verdict);
		assertTrue(warnings.isEmpty(), "switched off is not a warning");

		VariantFormRanking.Inputs noDiann = new VariantFormRanking.Inputs();
		noDiann.diann = new File(tmp, "no-diann.exe");
		VariantFormRanking.run(all, noDiann, tmp, s -> {
		}, warnings::add);
		assertEquals(VariantCandidate.Ms2Verdict.NOT_RANKED, open.ms2Verdict);
		assertEquals(1, warnings.size());
		assertFalse(new File(tmp, VariantFormRanking.WORK_SUBDIR).exists());
	}

	@Test
	@DisplayName("The verdict does not set the class; undecidable still wins")
	void classes() {
		VariantCandidate c = candidate(VariantCandidate.Engine.OPEN_SEARCH);
		c.inPangenomeAll = true;
		for (VariantCandidate.Ms2Verdict v : VariantCandidate.Ms2Verdict.values()) {
			c.ms2Verdict = v;
			c.classify();
			assertEquals(VariantCandidate.PangenomeClass.SUBSTITUTION_IN_PANGENOME, c.pangenomeClass, v.name());
		}
		c.ms2Verdict = VariantCandidate.Ms2Verdict.SUBSTITUTION;
		c.ms2Decidable = false;
		c.classify();
		assertEquals(VariantCandidate.PangenomeClass.UNDECIDABLE_BY_MS2, c.pangenomeClass);
	}
}
