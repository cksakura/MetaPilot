package uk.ac.quadram.metapilot.core.opensearch;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;

import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Test;

import uk.ac.quadram.metapilot.core.variant.OpenSearchVariantCandidateSource;
import uk.ac.quadram.metapilot.core.variant.VariantCandidate;

class OpenSearchVariantCandidateSourceTest {

	private static final double TOL = 0.005;

	private static OpenSearchPsm psm(String peptide, int scan, int peakId, boolean target) {
		OpenSearchPsm psm = new OpenSearchPsm();
		psm.setFileName("run1");
		psm.setScan(scan);
		psm.setHitRank(1);
		psm.setPeptide(peptide);
		psm.setPeakId(peakId);
		psm.setTarget(target);
		psm.setLocalizedSequence(peptide);
		return psm;
	}

	private static SubstitutionReassigner.Proposal proposal(OpenSearchPsm psm, int position, char from, char to) {
		double delta = IsobaricModificationTable.substitutionMass(from, to, SubstitutionMassTable.FIXED_CYSTEINE_MOD);
		StringBuilder candidate = new StringBuilder(psm.getPeptide());
		candidate.setCharAt(position, to);
		return new SubstitutionReassigner.Proposal(psm, candidate.toString(),
				new SubstitutionMassTable.Substitution(from, to, delta), position, true);
	}

	@Test
	@DisplayName("Shaped and rejected-peak proposals on validated targets become candidates; others are dropped")
	void keepsValidatedTargetsOnly() {
		OpenSearchPsm shaped = psm("PEPGTIDEK", 10, 3, true);
		shaped.setSubstitutionVerdict("substitution");
		shaped.setCharge(2);
		shaped.setPrecursorNeutralMass(1000.0);
		shaped.setQValue(0.002);
		OpenSearchPsm rejected = psm("SAMPLEK", 11, 7, true);
		OpenSearchPsm notValidated = psm("PEPGTIDEK", 12, 3, true);
		OpenSearchPsm decoy = psm("KEDITGPEP", 13, 3, false);

		List<SubstitutionReassigner.Proposal> proposals = new ArrayList<>(Arrays.asList(
				proposal(shaped, 3, 'G', 'A'), proposal(notValidated, 3, 'G', 'A'), proposal(decoy, 5, 'G', 'A')));
		List<SubstitutionReassigner.Proposal> rejectedProposals = new ArrayList<>(
				Arrays.asList(proposal(rejected, 0, 'S', 'T')));
		SubstitutionReassigner.Result result = new SubstitutionReassigner.Result(proposals, new ArrayList<>(),
				new HashSet<>(Arrays.asList(3)), rejectedProposals);

		OpenSearchVariantCandidateSource.Result out = OpenSearchVariantCandidateSource.collect(result,
				Arrays.asList(shaped, rejected, decoy), TOL);

		assertEquals(2, out.candidates.size());
		assertEquals(1, out.shapedCandidates);
		assertEquals(1, out.rejectedPeakCandidates);
		assertEquals(new HashSet<>(Arrays.asList("PEPGTIDEK", "SAMPLEK")), new HashSet<>(out.basePeptides));

		VariantCandidate a = out.candidates.get(0);
		assertEquals(VariantCandidate.Engine.OPEN_SEARCH, a.engine);
		assertEquals("PEPATIDEK", a.variantPeptide);
		assertEquals(4, a.position);
		assertEquals('G', a.fromAa);
		assertEquals('A', a.toAa);
		assertEquals(3, a.peakId);
		assertTrue(a.peakSubstitutionShaped);
		assertTrue(a.ms2Decidable);
		assertEquals("substitution", a.frontEndEvidence);
		assertEquals("PEPGTIDEK", a.localizedSequence);
		assertEquals(Boolean.TRUE, a.reassignerInDb);
		assertEquals(1, a.psmId, "hit rank");
		assertEquals(0.002, a.qValue, 1e-12);
		assertEquals(1000.0 / 2 + 1.007276466812, a.precursorMz, 1e-9);

		VariantCandidate b = out.candidates.get(1);
		assertEquals("TAMPLEK", b.variantPeptide);
		assertEquals(1, b.position);
		assertFalse(b.peakSubstitutionShaped);
		assertFalse(b.ms2Decidable, "S->T is isobaric with methylation of S");
		assertEquals("peak_not_substitution_shaped", b.frontEndEvidence);
	}

	@Test
	@DisplayName("No reassigner result gives no candidates but still fills the control pool")
	void nullResult() {
		OpenSearchPsm a = psm("PEPTIDEK", 1, 0, true);
		OpenSearchVariantCandidateSource.Result out = OpenSearchVariantCandidateSource.collect(null,
				Arrays.asList(a), TOL);
		assertTrue(out.candidates.isEmpty());
		assertEquals(1, out.basePeptides.size());
	}
}
