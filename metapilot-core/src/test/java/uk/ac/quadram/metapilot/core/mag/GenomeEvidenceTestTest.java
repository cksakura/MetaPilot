package uk.ac.quadram.metapilot.core.mag;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertNotNull;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.Set;

import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Test;

import uk.ac.quadram.metapilot.core.mag.GenomeEvidenceTest.GenomeRow;
import uk.ac.quadram.metapilot.core.mag.GenomeEvidenceTest.PeptideClass;

class GenomeEvidenceTestTest {

	private static final String R1 = "r1";
	private static final String R2 = "r2";

	/**
	 * Two runs. r1 searched GA, GB, GC, GD; r2 searched GB, GC, GD, GE. GA and
	 * GB have 160 observed peptides each (r1). GE has 20 observed peptides (r2),
	 * each one substitution from a GA peptide, so GE is a relative absent from
	 * the run where GA was found. Peptide {@code x} is in a GA and a GE protein
	 * and observed in both runs: unique to GA in r1, to GE in r2. Decoy {@code y}
	 * is the same on the reversed side, in two rev_GA proteins and one rev_GE
	 * protein. Three decoys pass the peptide FDR ({@code y}, one in rev_GC, one
	 * in rev_GD); 100 low-scoring unexplained decoys spread evenly over the
	 * reversed genomes fill the pi window.
	 */
	private static final class Fixture {
		final Random rnd = new Random(7);
		final Set<String> used = new HashSet<>();
		final LinkedHashMap<String, String> proteins = new LinkedHashMap<>();
		final Map<String, Set<String>> runGenomes = new LinkedHashMap<>();
		final GenomeEvidenceTest.Hits hits = new GenomeEvidenceTest.Hits();
		final List<String> ga = peptides(160);
		final List<String> gb = peptides(160);
		final List<String> ge = new ArrayList<>();
		final String x = peptide();
		final String y = peptide();
		final double yScore;

		Fixture() {
			for (int i = 0; i < 20; i++) {
				ge.add(variant(ga.get(i)));
			}
			runGenomes.put(R1, new HashSet<>(Arrays.asList("GA", "GB", "GC", "GD")));
			runGenomes.put(R2, new HashSet<>(Arrays.asList("GB", "GC", "GD", "GE")));

			proteins("GA", ga);
			proteins.put("GA_00033", x);
			proteins("GB", gb);
			proteins("GC", peptides(10));
			proteins("GD", peptides(10));
			proteins("GE", ge);
			proteins.put("GE_00005", x);

			String[] genomes = { "GA", "GB", "GC", "GD", "GE" };
			List<List<String>> low = new ArrayList<>();
			for (String g : genomes) {
				List<String> d = peptides(20);
				low.add(d);
				proteins("rev_" + g, d);
			}
			String hc = peptide();
			String hd = peptide();
			proteins.put("rev_GA_00001", proteins.get("rev_GA_00001") + y);
			proteins.put("rev_GA_00002", proteins.get("rev_GA_00002") + y);
			proteins.put("rev_GE_00001", proteins.get("rev_GE_00001") + y);
			proteins.put("rev_GC_00005", hc);
			proteins.put("rev_GD_00005", hd);

			List<String[]> targets = new ArrayList<>();
			for (String p : ga) {
				targets.add(new String[] { R1, p });
			}
			for (String p : gb) {
				targets.add(new String[] { R1, p });
			}
			for (String p : ge) {
				targets.add(new String[] { R2, p });
			}
			targets.add(new String[] { R1, x });
			targets.add(new String[] { R2, x });
			for (int i = 0; i < targets.size(); i++) {
				hits.add(targets.get(i)[0], targets.get(i)[1], true, PeptideClass.UNMODIFIED, 1000 - i);
			}
			hits.add(R1, hc, false, PeptideClass.UNMODIFIED, 1000 - 110.5);
			hits.add(R1, hd, false, PeptideClass.UNMODIFIED, 1000 - 220.5);
			yScore = 1000 - 320.5;
			hits.add(R1, y, false, PeptideClass.UNMODIFIED, yScore);
			hits.add(R2, y, false, PeptideClass.UNMODIFIED, yScore);
			for (int j = 0; j < 20; j++) {
				for (int k = 0; k < genomes.length; k++) {
					String run = "GE".equals(genomes[k]) ? R2 : R1;
					hits.add(run, low.get(k).get(j), false, PeptideClass.UNEXPLAINED, 100 - (j * 5 + k) * 0.1);
				}
			}
		}

		String peptide() {
			while (true) {
				StringBuilder sb = new StringBuilder();
				for (int i = 0; i < 9; i++) {
					sb.append("ACDEFGHLNPQSTVWY".charAt(rnd.nextInt(16)));
				}
				String p = sb.append('K').toString();
				if (used.add(p)) {
					return p;
				}
			}
		}

		List<String> peptides(int n) {
			List<String> out = new ArrayList<>();
			for (int i = 0; i < n; i++) {
				out.add(peptide());
			}
			return out;
		}

		String variant(String p) {
			for (char c : "ACDEFGHLNPQSTVWY".toCharArray()) {
				String v = p.substring(0, 3) + c + p.substring(4);
				if (c != p.charAt(3) && used.add(v)) {
					return v;
				}
			}
			throw new IllegalStateException(p);
		}

		/** Five peptides per protein, accessions {@code <prefix>_00001}, ... */
		void proteins(String prefix, List<String> peptides) {
			for (int i = 0; i * 5 < peptides.size(); i++) {
				proteins.put(String.format("%s_%05d", prefix, i + 1),
						String.join("", peptides.subList(i * 5, Math.min(peptides.size(), i * 5 + 5))));
			}
		}

		GenomeEvidenceTest.Result run(GenomeEvidenceTest.Settings settings) throws IOException {
			return GenomeEvidenceTest.run(hits, runGenomes, c -> proteins.forEach(c::accept), settings);
		}
	}

	private static Map<String, GenomeRow> byGenome(GenomeEvidenceTest.Result result) {
		Map<String, GenomeRow> out = new HashMap<>();
		for (GenomeRow row : result.rows) {
			out.put(row.genome, row);
		}
		return out;
	}

	@Test
	@DisplayName("Present genomes accepted, reversed genomes not")
	void acceptsPresentGenomesAndNoReversedGenome() throws IOException {
		GenomeEvidenceTest.Result result = new Fixture().run(null);
		Map<String, GenomeRow> rows = byGenome(result);

		assertTrue(rows.get("GA").accepted);
		assertTrue(rows.get("GB").accepted);
		assertFalse(rows.get("GC").accepted);
		assertFalse(rows.get("GD").accepted);
		assertEquals(160, rows.get("GA").keptPeptides);
		assertEquals(2, rows.get("GB").runsSearched);
		assertEquals(1, rows.get("GA").runsSearched);

		assertEquals(0, result.calibration.decoyGenomesAccepted);
		for (GenomeRow row : result.rows) {
			if (row.decoy) {
				assertFalse(row.accepted, row.genome);
				assertTrue(row.genome.startsWith("rev_"), row.genome);
			}
		}
		assertEquals(1, rows.get("rev_GC").keptPeptides);
		assertEquals(10, result.rows.size());
		assertEquals(2, result.calibration.decoyEvidencePeptides);
		assertEquals(60, result.calibration.windowSize);

		double[] stats = result.genomeStats.get("GA");
		assertEquals(160, stats[0], 0);
		assertTrue(stats[2] <= 0.01);
		assertEquals(1.0 - stats[2], result.genomeScoreMap.get("GA"), 1e-12);
		assertEquals(Set.of("GA", "GB"), DdaMagProteinInferenceService.confidentGenomes(result.genomeScoreMap,
				result.genomeStats, DdaMagProteinInferenceService.DEFAULT_GENOME_Q_VALUE_THRESHOLD));
	}

	@Test
	@DisplayName("A one-substitution relative absent from the run's database is excluded after its neighbour")
	void excludesNearNeighbourOfAcceptedGenome() throws IOException {
		Fixture fixture = new Fixture();
		GenomeRow ge = byGenome(fixture.run(null)).get("GE");
		assertEquals(20, ge.uniquePeptides);
		assertEquals(20, ge.excludedNearNeighbour);
		assertEquals(0, ge.keptPeptides);
		assertFalse(ge.accepted);

		GenomeEvidenceTest.Settings exact = new GenomeEvidenceTest.Settings();
		exact.maxSubstitutions = 0;
		GenomeRow geExact = byGenome(fixture.run(exact)).get("GE");
		assertEquals(20, geExact.keptPeptides);
		assertTrue(geExact.accepted);
	}

	@Test
	@DisplayName("A genome without decoys in the pi window still has a nonzero expected false count")
	void genomeWithoutWindowDecoysNeedsMoreThanOnePeptide() throws IOException {
		Fixture fixture = new Fixture();
		String single = fixture.peptide();
		fixture.proteins.put("GF_00001", single);
		fixture.runGenomes.get(R1).add("GF");
		fixture.hits.add(R1, single, true, PeptideClass.UNMODIFIED, 1000.5);

		GenomeEvidenceTest.Result result = fixture.run(null);
		GenomeRow gf = byGenome(result).get("GF");
		assertEquals(1, gf.keptPeptides);
		assertEquals(2.0 / (60 + 6), gf.expectedFalse, 1e-12);
		assertFalse(gf.accepted);
		assertTrue(byGenome(result).get("GA").accepted);
	}

	@Test
	@DisplayName("A peptide unique to different genomes in different runs is dropped, targets and decoys alike")
	void dropsPeptidesUniqueToDifferentGenomesInDifferentRuns() throws IOException {
		GenomeEvidenceTest.Result result = new Fixture().run(null);
		String inconsistent = "unique_to_different_genomes_in_different_runs";
		assertEquals(1, result.calibration.targetUniqueness.get(inconsistent));
		assertEquals(340, result.calibration.targetUniqueness.get("unique"));
		assertEquals(1, result.calibration.decoyUniqueness.get(inconsistent));
		assertEquals(2, result.calibration.decoyUniqueness.get("unique"));
		assertEquals(160, byGenome(result).get("GA").uniquePeptides);
	}

	@Test
	@DisplayName("Picked protein FDR: a decoy peptide shared by reversed proteins goes to one of them")
	void sharedDecoyPeptideScoresOneReversedProtein() throws IOException {
		Fixture fixture = new Fixture();
		GenomeEvidenceTest.Result result = fixture.run(null);
		double[] first = result.proteinStats.get("GA_00001");
		double[] second = result.proteinStats.get("GA_00002");
		assertNotNull(first);
		assertNotNull(second);
		assertEquals(1000, first[0], 0);
		assertEquals(fixture.yScore, first[1], 0);
		assertTrue(Double.isNaN(second[1]));
		assertEquals(0.0, first[2], 0);
		assertFalse(result.proteinStats.containsKey("GE_00001"));
		assertEquals(1, result.calibration.decoyProteinsScored);
	}

	@Test
	@DisplayName("Closed-stage decoys stay out of the pi window; the open-stage chi-square is reported")
	void closedDecoysAreNotWindowDecoys() throws IOException {
		Fixture fixture = new Fixture();
		List<String> closedDecoys = fixture.peptides(10);
		fixture.proteins.put("rev_GB_00099", String.join("", closedDecoys));
		for (String p : closedDecoys) {
			fixture.hits.add(R1, p, false, PeptideClass.CLOSED, 500);
		}
		fixture.hits.add(R1, fixture.ga.get(0), true, PeptideClass.CLOSED, 1000);

		GenomeEvidenceTest.Result result = fixture.run(null);

		assertEquals(102, result.calibration.windowAvailable);
		assertEquals(60, result.calibration.windowSize);
		assertEquals(2, result.calibration.decoyEvidencePeptides);
		assertEquals(2, result.calibration.openDecoyEvidencePeptides);
		assertEquals(result.calibration.chiSquare, result.calibration.chiSquareOpen, 1e-12);
		assertTrue(byGenome(result).get("GA").accepted);
	}

	@Test
	@DisplayName("Picked protein FDR: target against its reversed protein, ties to the decoy")
	void pickedProteinFdr() {
		Map<String, Double> targets = new HashMap<>();
		targets.put("P1", 10.0);
		targets.put("P2", 8.0);
		targets.put("P3", 5.0);
		Map<String, Double> decoys = new HashMap<>();
		decoys.put("rev_P2", 9.0);
		decoys.put("rev_P3", 5.0);
		decoys.put("rev_P4", 7.0);
		HashMap<String, double[]> out = GenomeEvidenceTest.pickedProteinFdr(targets, decoys, "rev_");
		assertEquals(3, out.size());
		assertEquals(0.0, out.get("P1")[2], 0);
		assertTrue(Double.isNaN(out.get("P1")[1]));
		assertEquals(9.0, out.get("P2")[1], 0);
		assertEquals(1.0, out.get("P2")[2], 0);
		assertEquals(1.0, out.get("P3")[2], 0);
	}

	@Test
	void logPoissonSurvival() {
		assertEquals(0.0, GenomeEvidenceTest.logPoissonSf(0, 3.0), 0);
		assertEquals(Math.log(1 - Math.exp(-0.4)), GenomeEvidenceTest.logPoissonSf(1, 0.4), 1e-12);
		assertEquals(Math.log(1 - 5 * Math.exp(-2)), GenomeEvidenceTest.logPoissonSf(3, 2.0), 1e-12);
		double big = GenomeEvidenceTest.logPoissonSf(500, 1.0);
		assertTrue(Double.isFinite(big));
		assertTrue(big < -2000);
		assertTrue(GenomeEvidenceTest.logPoissonSf(501, 1.0) < big);
	}

	@Test
	void withinSubstitutions() {
		assertTrue(GenomeEvidenceTest.withinSubstitutions("xxPEPTLDEK", 2, "PEPTLDEK", 0));
		assertTrue(GenomeEvidenceTest.withinSubstitutions("PEPTLDEK", 0, "PEPTNDEK", 1));
		assertFalse(GenomeEvidenceTest.withinSubstitutions("PEPTLDEK", 0, "PAPTNDEK", 1));
		assertTrue(GenomeEvidenceTest.withinSubstitutions("PEPTLDEK", 0, "PAPTNDEK", 2));
	}

	@Test
	@DisplayName("RESET decoy weight counts each passed decoy twice in the peptide walk and the expected false count")
	void decoyWeightCountsEachDecoyTwice() throws IOException {
		GenomeEvidenceTest.Result plain = new Fixture().run(null);
		GenomeEvidenceTest.Settings reset = new GenomeEvidenceTest.Settings();
		reset.decoyWeight = 2.0;
		GenomeEvidenceTest.Result weighted = new Fixture().run(reset);
		int u = PeptideClass.UNMODIFIED.ordinal();
		// Plain: all 341 target units and the 3 decoys pass (3/341 < 1%). Weighted: 2 x 1/221 < 1% just before
		// the second decoy, 2 x 2/221 > 1% after it, and the walk never comes back under 1%.
		assertEquals(341, plain.calibration.acceptedUnits[u][0]);
		assertEquals(3, plain.calibration.acceptedUnits[u][1]);
		assertEquals(221, weighted.calibration.acceptedUnits[u][0]);
		assertEquals(1, weighted.calibration.acceptedUnits[u][1]);
		// Same F and window: the expected false count of a target genome doubles.
		GenomeEvidenceTest.Settings sameWalk = new GenomeEvidenceTest.Settings();
		sameWalk.decoyWeight = 2.0;
		sameWalk.peptideFdr = 0.05;
		GenomeEvidenceTest.Settings plainWide = new GenomeEvidenceTest.Settings();
		plainWide.peptideFdr = 0.05;
		Map<String, GenomeRow> a = byGenome(new Fixture().run(plainWide));
		Map<String, GenomeRow> b = byGenome(new Fixture().run(sameWalk));
		assertEquals(2 * a.get("GB").expectedFalse, b.get("GB").expectedFalse, 1e-9);
	}

	@Test
	@DisplayName("Protein filter keeps proteins with two or more peptides and walks one-peptide proteins to the FDR")
	void proteinFilterWalksOnePeptideProteins() {
		Map<String, Double> targetScores = new HashMap<>();
		Map<String, Integer> targetPeptides = new HashMap<>();
		for (int i = 0; i < 200; i++) {
			targetScores.put("T" + i, 300.0 - i);
			targetPeptides.put("T" + i, 1);
		}
		for (int i = 0; i < 3; i++) {
			targetScores.put("M" + i, 1.0);
			targetPeptides.put("M" + i, 2);
		}
		Map<String, Double> decoyScores = new HashMap<>();
		Map<String, Integer> decoyPeptides = new HashMap<>();
		decoyScores.put("rev_D0", 150.5);
		decoyPeptides.put("rev_D0", 1);
		decoyScores.put("rev_D1", 0.0);
		decoyPeptides.put("rev_D1", 1);
		decoyScores.put("rev_D2", 500.0);
		decoyPeptides.put("rev_D2", 3);
		// Weight 1 at 1.4%: (1 + 1) / 200 = 1% above the last decoy, so every target passes; that decoy is at 3 / 200.
		GenomeEvidenceTest.Calibration plain = new GenomeEvidenceTest.Calibration();
		Set<String> all = GenomeEvidenceTest.proteinFilter(targetScores, targetPeptides, decoyScores, decoyPeptides,
				0.014, 1.0, plain);
		assertEquals(203, all.size());
		assertEquals(3, plain.proteinsTwoOrMore[0]);
		assertEquals(1, plain.proteinsTwoOrMore[1], "multi-peptide decoys are counted, never walked");
		assertEquals(200, plain.proteinsOnePeptide[0]);
		assertEquals(2, plain.proteinsOnePeptide[1]);
		assertEquals(1, plain.proteinsOnePeptidePassing[1]);
		// Weight 2: 2 x 1 / 150 = 1.3% above the first decoy; 2 x 2 / 200 = 2% below it.
		GenomeEvidenceTest.Calibration reset = new GenomeEvidenceTest.Calibration();
		Set<String> kept = GenomeEvidenceTest.proteinFilter(targetScores, targetPeptides, decoyScores, decoyPeptides,
				0.014, 2.0, reset);
		assertEquals(153, kept.size());
		assertTrue(kept.contains("T149") && !kept.contains("T150") && kept.contains("M2"));
		assertEquals(0, reset.proteinsOnePeptidePassing[1]);
	}

	@Test
	@DisplayName("Host proteins are no candidate genome, and a peptide they contain is no genome evidence")
	void hostPeptidesAreNotGenomeEvidence() throws IOException {
		Fixture f = new Fixture();
		// Ten of GB's peptides also occur in a host protein; one host protein shares nothing.
		f.proteins.put("sp|P00001|HOST1_HUMAN", String.join("", f.gb.subList(0, 10)));
		f.proteins.put("sp|P00002|HOST2_HUMAN", String.join("", f.peptides(5)));
		GenomeEvidenceTest.Settings settings = new GenomeEvidenceTest.Settings();
		settings.genomePrefix = "G";
		GenomeEvidenceTest.Result result = f.run(settings);
		Map<String, GenomeRow> rows = byGenome(result);

		assertEquals(150, rows.get("GB").uniquePeptides);
		assertEquals(10, result.calibration.targetUniqueness.get("host").intValue());
		assertTrue(rows.get("GB").accepted);
		for (GenomeRow row : result.rows) {
			assertTrue(row.genome.startsWith("G") || row.genome.startsWith("rev_G"), row.genome);
		}

		// Without the prefix the host is ignored as before, and the shared peptides count for GB.
		assertEquals(160, byGenome(f.run(null)).get("GB").uniquePeptides);
	}
}
