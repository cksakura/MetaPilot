package uk.ac.quadram.metapilot.core.variant;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.io.File;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

import org.junit.jupiter.api.Assumptions;
import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Test;

/**
 * Agreement with the reference implementation {@code variant_pangenome_check.py}
 * shipped with the plan, on real data. Skipped unless pointed at its output:
 *
 * <pre>
 *   python variant_pangenome_check.py pFind-Filtered.spectra uhgp-100.faa.gz out
 *   mvn -pl metapilot-core test -Dtest=PythonOracleComparisonTest \
 *       -Dpangenome.oracle.peptides=out_peptides.tsv -Dpangenome.oracle.fasta=uhgp-100.faa.gz \
 *       [-Dpangenome.oracle.spectra=pFind-Filtered.spectra]
 * </pre>
 *
 * The Python script has no tryptic-context rule and no species scope, so the
 * Java lookup runs with both off. Peptides containing a letter outside the 20
 * residues are left out of the comparison. Two checks:
 * <ol>
 * <li>Every peptide the script looked up (variants, controls and bases) gets the
 * same found / not-found answer from the Java lookup over the same FASTA.</li>
 * <li>With the spectra file given: every single-substitution candidate the Java
 * pFind front end produces from a plain token is among the script's variant
 * peptides (the script does not parse pFind's alias form). The
 * reverse does not hold by design - the Java reader drops peptides shorter
 * than 7, sequences with X, PSMs with any decoy protein and K/R
 * substitutions, which the script keeps.</li>
 * </ol>
 */
class PythonOracleComparisonTest {

	@Test
	@DisplayName("Java lookup reproduces variant_pangenome_check.py in_pangenome per peptide")
	void agreesWithPythonOracle() throws IOException {
		String peptidesPath = System.getProperty("pangenome.oracle.peptides");
		String fastaPath = System.getProperty("pangenome.oracle.fasta");
		Assumptions.assumeTrue(peptidesPath != null && fastaPath != null,
				"set -Dpangenome.oracle.peptides and -Dpangenome.oracle.fasta to run");

		Map<String, Boolean> oracle = new LinkedHashMap<>();
		Map<String, String> classes = new LinkedHashMap<>();
		List<String> lines = Files.readAllLines(new File(peptidesPath).toPath(), StandardCharsets.UTF_8);
		for (String line : lines.subList(1, lines.size())) {
			String[] f = line.split("\t", -1);
			String peptide = f[0].replace('J', 'L');
			// The script would match an X in a query against an X in a protein; the Java lookup never
			// matches across a non-residue. Such peptides are not substitution candidates either way.
			if (!peptide.matches("[ACDEFGHKLMNPQRSTVWY]+")) {
				continue;
			}
			oracle.merge(peptide, "1".equals(f[2]), (a, b) -> a || b);
			classes.merge(peptide, f[1], (a, b) -> a + "," + b);
		}
		List<String> queries = new ArrayList<>(oracle.keySet());
		PangenomeVariantLookup.QueryHit[] hits = PangenomeVariantLookup.scan(queries,
				List.of(new File(fastaPath)), null, false, Math.max(1, Runtime.getRuntime().availableProcessors() - 1),
				"", System.out::println, "oracle FASTA");
		int disagreements = 0;
		StringBuilder examples = new StringBuilder();
		for (int i = 0; i < queries.size(); i++) {
			boolean expected = oracle.get(queries.get(i));
			if (hits[i].found() != expected) {
				disagreements++;
				if (disagreements <= 20) {
					examples.append('\n').append(queries.get(i)).append(" (").append(classes.get(queries.get(i)))
							.append("): python=").append(expected).append(" java=").append(hits[i].found());
				}
			}
		}
		assertEquals(0, disagreements, disagreements + " of " + queries.size() + " peptides disagree:" + examples);

		String spectraPath = System.getProperty("pangenome.oracle.spectra");
		if (spectraPath == null) {
			return;
		}
		PFindVariantCandidateSource.Result java = PFindVariantCandidateSource.read(List.of(new File(spectraPath)), "",
				0.005, System.out::println);
		int missing = 0;
		StringBuilder missingExamples = new StringBuilder();
		for (VariantCandidate c : java.candidates) {
			// The script's token regex skips pFind's alias form, Methyl[S](Ser->Thr[S]).
			boolean alias = c.isobaricModifications.stream().anyMatch(s -> s.startsWith("pFind alias"));
			if (alias || differences(c.basePeptide, c.variantPeptide) != 1) {
				continue;
			}
			String v = PeptideAutomaton.normalise(c.variantPeptide);
			String cl = classes.get(v);
			if (cl == null || !cl.contains("real_")) {
				missing++;
				if (missing <= 20) {
					missingExamples.append('\n').append(c.basePeptide).append(' ').append(c.frontEndEvidence);
				}
			}
		}
		assertTrue(missing == 0, missing + " Java candidates are not among the script's variants:" + missingExamples);
	}

	private static int differences(String a, String b) {
		String x = PeptideAutomaton.normalise(a);
		String y = PeptideAutomaton.normalise(b);
		int n = 0;
		for (int i = 0; i < x.length(); i++) {
			if (x.charAt(i) != y.charAt(i)) {
				n++;
			}
		}
		return n;
	}
}
