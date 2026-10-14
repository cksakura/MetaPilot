package uk.ac.quadram.metapilot.core.opensearch;

import static org.junit.jupiter.api.Assertions.assertEquals;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;

import org.junit.jupiter.api.Test;

class UnmodifiedPsmSelectorTest {

	@Test
	void reportsTargetDecoyCountsAndAchievedFdr() {
		List<OpenSearchPsm> psms = new ArrayList<>();
		Random random = new Random(23L);
		for (int i = 0; i < 100; i++) {
			psms.add(psm(true, random.nextGaussian() * 0.0025, 1e-10));
		}
		psms.add(psm(false, random.nextGaussian() * 0.0025, 1e-5));

		UnmodifiedPsmSelector.Result result = UnmodifiedPsmSelector.select(psms,
				new DeltaMassPeakFitter.Config(), 0.01, null);

		assertEquals(100, result.targetsBeforeFdr);
		assertEquals(1, result.decoysBeforeFdr);
		assertEquals(100, result.targetsAfterFdr);
		assertEquals(1, result.decoysAfterFdr);
		assertEquals(0.01, result.estimatedFdr, 1e-12);
		assertEquals(101, result.unmodified.size());
		assertEquals(0.01, result.unmodified.get(100).getQValue(), 1e-12);
	}

	private static OpenSearchPsm psm(boolean target, double massDiff, double expect) {
		OpenSearchPsm psm = new OpenSearchPsm();
		psm.setTarget(target);
		psm.setMassDiff(massDiff);
		psm.setExpect(expect);
		return psm;
	}
}
