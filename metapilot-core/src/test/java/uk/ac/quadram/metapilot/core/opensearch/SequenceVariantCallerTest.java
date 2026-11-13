package uk.ac.quadram.metapilot.core.opensearch;

import static org.junit.jupiter.api.Assertions.assertEquals;

import java.util.ArrayList;
import java.util.List;

import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Test;

class SequenceVariantCallerTest {

	private static OpenSearchPsm psm(String file, int scan, String peptide, double delta) {
		OpenSearchPsm p = new OpenSearchPsm();
		p.setFileName(file);
		p.setScan(scan);
		p.setPeptide(peptide);
		p.setTarget(true);
		p.setHitRank(1);
		p.setCharge(2);
		p.setMassDiff(delta);
		p.setNormalizedMassDiff(delta);
		p.setPeakId(-1);
		p.setProtein("MGYG000000001_00001");
		return p;
	}

	/** A->T (one nucleotide) on SAMPLEK, set aside in three samples. */
	private static List<OpenSearchPsm> evidence() {
		List<OpenSearchPsm> out = new ArrayList<>();
		for (int i = 1; i <= 3; i++) {
			OpenSearchPsm p = psm("run" + i, 10 + i, "SAMPLEK", 30.0106);
			p.setSubstitutionVerdict(SubstitutionReassigner.Verdict.SUBSTITUTION_NOT_IN_DB.name());
			p.setSubstitutionCandidate("STMPLEK");
			out.add(p);
		}
		return out;
	}

	private static SequenceVariantCaller.VariantClass classOf(SequenceVariantCaller.Result r) {
		assertEquals(1, r.variants.size());
		return r.variants.get(0).variantClass;
	}

	@Test
	@DisplayName("The reference is 'also observed' when any identification shows it unmodified, the closed "
			+ "search's included; the variant evidence alone cannot say so")
	void referenceObservationsComeFromTheWholeResult() {
		SequenceVariantCaller.Config cfg = new SequenceVariantCaller.Config();
		// Without reference observations given, only the evidence list is read: nothing unmodified in it.
		assertEquals(SequenceVariantCaller.VariantClass.SEQUENCE_VARIANT,
				classOf(SequenceVariantCaller.call(evidence(), List.of(), cfg, null)));
		// An exact-mass identification of the reference peptide, as the closed search reports it.
		OpenSearchPsm closed = psm("run2", 99, "SAMPLEK", 0.0);
		closed.setClosedStage(true);
		SequenceVariantCaller.Result withReference = SequenceVariantCaller.call(evidence(), List.of(closed), List.of(),
				cfg, null);
		assertEquals(SequenceVariantCaller.VariantClass.SEQUENCE_VARIANT_WITH_REFERENCE, classOf(withReference));
		assertEquals(1, withReference.variants.get(0).coObservedSamples());
		// Given but empty: still no reference, and the evidence list is not searched for one.
		assertEquals(SequenceVariantCaller.VariantClass.SEQUENCE_VARIANT,
				classOf(SequenceVariantCaller.call(evidence(), new ArrayList<>(), List.of(), cfg, null)));
		// Another peptide does not count.
		assertEquals(SequenceVariantCaller.VariantClass.SEQUENCE_VARIANT, classOf(SequenceVariantCaller.call(evidence(),
				List.of(psm("run1", 5, "SAMPLER", 0.0)), List.of(), cfg, null)));
	}
}
