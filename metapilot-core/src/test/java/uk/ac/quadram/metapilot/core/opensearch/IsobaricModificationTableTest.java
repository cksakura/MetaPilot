package uk.ac.quadram.metapilot.core.opensearch;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.util.List;
import java.util.Set;

import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Test;

class IsobaricModificationTableTest {

	private static final double TOL = 0.005;

	private static Set<String> same(char from, char to) {
		double delta = IsobaricModificationTable.substitutionMass(from, to, 0.0);
		return IsobaricModificationTable.sameResidue(from, delta, false, TOL);
	}

	@Test
	@DisplayName("N->D and Q->E are deamidation of the same residue")
	void deamidation() {
		assertTrue(same('N', 'D').stream().anyMatch(s -> s.startsWith("Deamidation")), same('N', 'D').toString());
		assertTrue(same('Q', 'E').stream().anyMatch(s -> s.startsWith("Deamidation")), same('Q', 'E').toString());
	}

	@Test
	@DisplayName("D->E, S->T and N->Q are Unimod:34 methylation of the same residue")
	void methylation() {
		assertTrue(same('D', 'E').contains("Methyl (D) [Unimod:34]"));
		assertTrue(same('S', 'T').contains("Methyl (S) [Unimod:34]"));
		assertTrue(same('N', 'Q').contains("Methyl (N) [Unimod:34]"));
	}

	@Test
	@DisplayName("G->A and V->L are decidable: no methylation is listed on G or V")
	void glycineAndValineDecidable() {
		assertTrue(same('G', 'A').isEmpty(), same('G', 'A').toString());
		assertTrue(same('V', 'L').isEmpty(), same('V', 'L').toString());
	}

	@Test
	@DisplayName("A->S is decidable: oxidation sits on M, not on A")
	void oxidationOnlyOnItsOwnResidue() {
		assertTrue(same('A', 'S').isEmpty(), same('A', 'S').toString());
		double oxidation = IsobaricModificationTable.substitutionMass('A', 'S', 0.0);
		assertTrue(IsobaricModificationTable.sameResidue('M', oxidation, false, TOL).contains("Oxidation (M)"));
	}

	@Test
	@DisplayName("Cysteine mass follows the engine's fixed modification")
	void cysteineConvention() {
		double unmodified = IsobaricModificationTable.substitutionMass('S', 'C', 0.0);
		double carbamidomethyl = IsobaricModificationTable.substitutionMass('S', 'C',
				SubstitutionMassTable.FIXED_CYSTEINE_MOD);
		assertEquals(SubstitutionMassTable.FIXED_CYSTEINE_MOD, carbamidomethyl - unmodified, 1e-9);
	}

	@Test
	@DisplayName("The logged table lists the flagged pairs and leaves out K/R and I<->L")
	void describe() {
		List<String> table = IsobaricModificationTable.describe(TOL, 0.0);
		assertTrue(table.stream().anyMatch(s -> s.startsWith("N->D ")));
		assertTrue(table.stream().anyMatch(s -> s.startsWith("D->E ")));
		assertFalse(table.stream().anyMatch(s -> s.startsWith("G->A ")));
		assertFalse(table.stream().anyMatch(s -> s.startsWith("K->") || s.contains("->K ") || s.startsWith("R->")
				|| s.contains("->R ")));
		assertFalse(table.stream().anyMatch(s -> s.startsWith("I->L ") || s.startsWith("L->I ")));
	}
}
