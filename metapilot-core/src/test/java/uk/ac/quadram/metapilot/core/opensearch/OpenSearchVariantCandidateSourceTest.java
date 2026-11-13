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

	private static ShiftLocaliser.Result localised(int site, String... explanations) {
		return new ShiftLocaliser.Result(7.5, site, 0.95, ShiftLocaliser.State.LOCALISED_EXPLAINED,
				Arrays.asList(explanations));
	}

	@Test
	@DisplayName("Localised substitutions on validated targets become candidates at their site; a matching peak "
			+ "proposal gains the localisation instead of a second candidate")
	void localisedCandidates() {
		OpenSearchPsm offPeak = psm("PEPGTIDEK", 20, -1, true);
		offPeak.setShiftLocalisation(localised(3, "subst G->A", "unimod Methyl"));
		OpenSearchPsm onPeak = psm("SAMGLEK", 21, 5, true);
		onPeak.setSubstitutionVerdict("substitution");
		onPeak.setShiftLocalisation(localised(3, "subst G->A"));
		OpenSearchPsm ambiguous = psm("PEPGTIDEK", 22, -1, true);
		ambiguous.setShiftLocalisation(new ShiftLocaliser.Result(0.5, 3, 0.4, ShiftLocaliser.State.AMBIGUOUS,
				List.of()));
		OpenSearchPsm decoy = psm("KEDITGPEP", 23, -1, false);
		decoy.setShiftLocalisation(localised(5, "subst G->A"));
		OpenSearchPsm notValidated = psm("PEPGTIDEK", 24, -1, true);
		notValidated.setShiftLocalisation(localised(3, "subst G->A"));
		OpenSearchPsm modificationOnly = psm("PEPMTIDEK", 25, -1, true);
		modificationOnly.setShiftLocalisation(localised(3, "unimod Oxidation"));
		OpenSearchPsm acceptedShift = psm("PEPGTIDEK", 26, -1, true);
		acceptedShift.setShiftLocalisation(localised(3, "subst G->A"));
		acceptedShift.setIdentificationGroup("OPEN_ACCEPTED_SHIFT"); // an accepted-peak modification, not a variant
		for (OpenSearchPsm p : Arrays.asList(offPeak, onPeak, ambiguous, decoy, notValidated, modificationOnly)) {
			p.setIdentificationGroup("OPEN_LOCALISED");
		}
		offPeak.setIdentificationGroup("OPEN_SUBSTITUTION_SITE"); // reported, not identified, still a variant source
		onPeak.setIdentificationGroup("OPEN_SUBSTITUTION_SITE");

		// A peak proposal counts only where it is also the identification's call: the one for the
		// accepted-shift PSM, and one at another residue than the fragments' site, are dropped.
		SubstitutionReassigner.Result reassigner = new SubstitutionReassigner.Result(
				new ArrayList<>(Arrays.asList(proposal(onPeak, 3, 'G', 'A'), proposal(acceptedShift, 3, 'G', 'A'),
						proposal(offPeak, 0, 'P', 'S'))), new ArrayList<>(),
				new HashSet<>(Arrays.asList(5)), new ArrayList<>());
		OpenSearchVariantCandidateSource.Result out = OpenSearchVariantCandidateSource.collect(reassigner,
				Arrays.asList(offPeak, onPeak, ambiguous, decoy, modificationOnly, acceptedShift), TOL);

		assertEquals(2, out.candidates.size());
		assertEquals(1, out.shapedCandidates);
		assertEquals(1, out.localisedCandidates);
		VariantCandidate peak = out.candidates.get(0);
		assertEquals("SAMALEK", peak.variantPeptide);
		assertTrue(peak.peakSubstitutionShaped);
		assertEquals(0.95, peak.localisationSiteProbability, 1e-12, "the peak proposal gains the localisation");
		assertEquals(7.5, peak.localisationSupport, 1e-12);
		VariantCandidate loc = out.candidates.get(1);
		assertEquals("PEPATIDEK", loc.variantPeptide);
		assertEquals(4, loc.position, "1-based");
		assertEquals('G', loc.fromAa);
		assertEquals('A', loc.toAa);
		assertEquals(-1, loc.peakId);
		assertEquals(null, loc.peakSubstitutionShaped);
		assertEquals(null, loc.reassignerInDb);
		assertEquals(OpenSearchVariantCandidateSource.LOCALISED_CALL, loc.frontEndEvidence);
		assertEquals(0.95, loc.localisationSiteProbability, 1e-12);
		assertEquals(SubstitutionMassTable.residueMass('A') - SubstitutionMassTable.residueMass('G'), loc.deltaMass,
				1e-9);
		assertTrue(loc.ms2Decidable, "G->A has no same-residue modification of equal mass");

		// Localisation alone still gives candidates when the reassigner did not run.
		OpenSearchVariantCandidateSource.Result alone = OpenSearchVariantCandidateSource.collect(null,
				Arrays.asList(offPeak), TOL);
		assertEquals(1, alone.localisedCandidates);
		assertEquals("PEPATIDEK", alone.candidates.get(0).variantPeptide);
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
