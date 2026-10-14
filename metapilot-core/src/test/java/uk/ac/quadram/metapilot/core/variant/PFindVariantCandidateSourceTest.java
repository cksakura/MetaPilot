package uk.ac.quadram.metapilot.core.variant;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertNull;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.util.List;

import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Test;

class PFindVariantCandidateSourceTest {

	private static final double TOL = 0.005;

	private static List<VariantCandidate> candidates(String sequence, String modification,
			PFindVariantCandidateSource.Result stats) {
		return PFindVariantCandidateSource.fromPsm("raw1", 100, 2, sequence, modification, TOL, stats);
	}

	@Test
	@DisplayName("Entries split on the first comma only; non-substitution entries are ignored")
	void parseFirstComma() {
		List<PFindVariantCandidateSource.Token> t = PFindVariantCandidateSource
				.parse("3,Oxidation[M];5,Ser->Thr[S];0,Cation:Fe[III],x[AnyN-term];");
		assertEquals(1, t.size());
		assertEquals(5, t.get(0).position);
		assertEquals('S', t.get(0).from);
		assertEquals('T', t.get(0).to);
		assertNull(t.get(0).alias);
	}

	@Test
	@DisplayName("A plain substitution gives one decidable candidate with the variant applied")
	void plainSubstitution() {
		PFindVariantCandidateSource.Result stats = new PFindVariantCandidateSource.Result();
		List<VariantCandidate> c = candidates("PEPGIDEK", "4,Gly->Ala[G];", stats);
		assertEquals(1, c.size());
		VariantCandidate v = c.get(0);
		assertEquals("PEPAIDEK", v.variantPeptide);
		assertEquals(4, v.position);
		assertEquals('G', v.fromAa);
		assertEquals('A', v.toAa);
		assertEquals(14.01565, v.deltaMass, 1e-4);
		assertTrue(v.ms2Decidable);
		assertEquals(2, v.psmId);
		assertEquals(VariantCandidate.Engine.PFIND, v.engine);
		assertNull(v.peakSubstitutionShaped);
	}

	@Test
	@DisplayName("pFind's alias form marks the candidate undecidable and names the alias")
	void aliasUndecidable() {
		PFindVariantCandidateSource.Result stats = new PFindVariantCandidateSource.Result();
		List<VariantCandidate> c = candidates("PEPSIDEK", "4,Methyl[S](Ser->Thr[S]);", stats);
		assertEquals(1, c.size());
		assertFalse(c.get(0).ms2Decidable);
		assertTrue(c.get(0).isobaricModifications.contains("pFind alias Methyl[S]"),
				c.get(0).isobaricModifications.toString());
		assertEquals("PEPTIDEK", c.get(0).variantPeptide);
	}

	@Test
	@DisplayName("A same-residue modification of equal mass makes a plain call undecidable")
	void sameResidueUndecidable() {
		PFindVariantCandidateSource.Result stats = new PFindVariantCandidateSource.Result();
		List<VariantCandidate> c = candidates("PEPNIDEK", "4,Asn->Asp[N];", stats);
		assertEquals(1, c.size());
		assertFalse(c.get(0).ms2Decidable);
	}

	@Test
	@DisplayName("Xle matches either I or L in the sequence; I<->L sites are dropped")
	void xle() {
		PFindVariantCandidateSource.Result stats = new PFindVariantCandidateSource.Result();
		List<VariantCandidate> c = candidates("PEPIIDEK", "4,Xle->Val[I];5,Ile->Leu[I];", stats);
		assertEquals(1, c.size());
		assertEquals('I', c.get(0).fromAa);
		assertEquals('V', c.get(0).toAa);
		assertEquals("PEPVIDEK", c.get(0).variantPeptide);
	}

	@Test
	@DisplayName("A token naming the wrong residue skips the PSM")
	void residueMismatch() {
		PFindVariantCandidateSource.Result stats = new PFindVariantCandidateSource.Result();
		assertTrue(candidates("PEPGIDEK", "3,Gly->Ala[G];", stats).isEmpty());
		assertEquals(1, stats.skippedResidueMismatch);
		assertTrue(candidates("PEPGIDEK", "9,Gly->Ala[G];", stats).isEmpty());
		assertEquals(2, stats.skippedResidueMismatch);
	}

	@Test
	@DisplayName("A K or R substitution skips the whole PSM")
	void lysineArginine() {
		PFindVariantCandidateSource.Result stats = new PFindVariantCandidateSource.Result();
		assertTrue(candidates("PEPGIDEK", "4,Gly->Ala[G];8,Lys->Arg[K];", stats).isEmpty());
		assertEquals(1, stats.skippedKR);
	}

	@Test
	@DisplayName("Several substitutions give one candidate per site, all with every substitution applied")
	void multiSubstitution() {
		PFindVariantCandidateSource.Result stats = new PFindVariantCandidateSource.Result();
		List<VariantCandidate> c = candidates("PEPGIDEK", "1,Pro->Ala[P];4,Gly->Ser[G];", stats);
		assertEquals(2, c.size());
		assertEquals("AEPSIDEK", c.get(0).variantPeptide);
		assertEquals("AEPSIDEK", c.get(1).variantPeptide);
		assertEquals(1, c.get(0).position);
		assertEquals(4, c.get(1).position);
	}

	@Test
	@DisplayName("Cysteine is unmodified in a pFind open search")
	void cysteineUnmodified() {
		PFindVariantCandidateSource.Result stats = new PFindVariantCandidateSource.Result();
		List<VariantCandidate> c = candidates("PEPSIDEK", "4,Ser->Cys[S];", stats);
		assertEquals(103.00919 - 87.03203, c.get(0).deltaMass, 1e-4);
	}
}
