package uk.ac.quadram.metapilot.core.runner;

import uk.ac.quadram.metapilot.core.mag.GenomeEvidenceTest;
import uk.ac.quadram.metapilot.core.opensearch.OpenSearchLayout;
import uk.ac.quadram.metapilot.core.prodb.TargetDecoyFasta;
import uk.ac.quadram.metapilot.core.workflow.ProgressCallback;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;

/**
 * Feeds {@link GenomeEvidenceTest} from an open-search result: the spectrum
 * winners of {@code scored_psms.tsv} (one per spectrum, every class and score,
 * decoys included) and the database each run's stage 2 searched
 * ({@link OpenSearchLayout#refinedFasta}: the selected genomes' full proteomes
 * plus one {@code rev_} decoy per target).
 * <p>
 * Classes: unmodified when the mass difference is within
 * {@link #UNMODIFIED_TOLERANCE_DA} of 0, +1 or +2 isotope spacings; known shift
 * when the PSM sits on an accepted delta-mass peak; otherwise unexplained.
 * With the closed-then-open identification, closed-stage winners are the closed
 * class, and the open winner of a spectrum the closed stage accepted is in no
 * class (its decoys still join the {@code pi_g} window).
 */
public final class OpenSearchGenomeEvidence {

	static final double ISOTOPE_SPACING = 1.0033548;
	static final double UNMODIFIED_TOLERANCE_DA = 0.03;
	static final int MAX_ISOTOPE_ERROR = 2;

	private OpenSearchGenomeEvidence() {
	}

	/**
	 * @param resultRoot the engine search root holding {@code search_results/}
	 * @param engineMode the engine configuration whose {@code scored_psms.tsv} is read
	 * @param samples    sample base names, as in the PSM table's {@code file} column
	 */
	public static GenomeEvidenceTest.Result run(File resultRoot, String engineMode, List<String> samples,
			ProgressCallback cb) throws IOException {
		Map<String, File> databases = new LinkedHashMap<>();
		for (String sample : samples) {
			databases.put(sample, OpenSearchLayout.refinedFasta(resultRoot, sample));
		}
		GenomeEvidenceTest.Settings settings = new GenomeEvidenceTest.Settings();
		settings.decoyPrefix = TargetDecoyFasta.DECOY_TAG;
		return run(OpenSearchLayout.scoredPsms(resultRoot, engineMode), databases, settings, cb);
	}

	static GenomeEvidenceTest.Result run(File scoredPsms, Map<String, File> databases,
			GenomeEvidenceTest.Settings settings, ProgressCallback cb) throws IOException {
		if (cb == null) {
			cb = new ProgressCallback.NoOp();
		}
		cb.onMessage("Genome test: reading spectrum winners from " + scoredPsms.getName());
		GenomeEvidenceTest.Hits hits = readHits(scoredPsms);
		Map<String, Set<String>> runGenomes = new LinkedHashMap<>();
		for (String run : hits.runs()) {
			File fasta = databases.get(run);
			if (fasta == null || !fasta.isFile()) {
				throw new IOException("Genome test: no stage-2 database for run " + run + " ("
						+ (fasta == null ? "not a sample of this analysis" : fasta.getAbsolutePath() + " missing")
						+ "); peptide uniqueness is decided against each run's own database");
			}
			runGenomes.put(run, genomesOf(fasta, settings.decoyPrefix));
		}
		List<File> fastas = new ArrayList<>();
		Set<String> candidates = new HashSet<>();
		for (Map.Entry<String, Set<String>> e : runGenomes.entrySet()) {
			fastas.add(databases.get(e.getKey()));
			candidates.addAll(e.getValue());
		}
		cb.onMessage("Genome test: " + runGenomes.size() + " runs, " + candidates.size()
				+ " candidate genomes; digesting the runs' databases");
		GenomeEvidenceTest.Result result = GenomeEvidenceTest.run(hits, runGenomes, consumer -> forEachProtein(fastas,
				consumer), settings);
		GenomeEvidenceTest.Calibration c = result.calibration;
		cb.onMessage("Genome test: " + c.targetGenomesAccepted + " of " + c.candidateGenomes + " genomes accepted at q <= "
				+ c.genomeFdr + "; reversed genomes accepted " + c.decoyGenomesAccepted + "; F = "
				+ c.decoyEvidencePeptides + ", decoy window " + c.windowSize + " (chi-square p = " + c.chiSquareP + ")");
		return result;
	}

	static GenomeEvidenceTest.Hits readHits(File scoredPsms) throws IOException {
		GenomeEvidenceTest.Hits hits = new GenomeEvidenceTest.Hits();
		try (BufferedReader br = Files.newBufferedReader(scoredPsms.toPath(), StandardCharsets.UTF_8)) {
			String header = br.readLine();
			if (header == null) {
				throw new IOException(scoredPsms + " is empty");
			}
			Map<String, Integer> col = new HashMap<>();
			String[] names = header.split("\t", -1);
			for (int i = 0; i < names.length; i++) {
				col.put(names[i], i);
			}
			int file = column(col, "file", scoredPsms);
			int peptide = column(col, "peptide", scoredPsms);
			int target = column(col, "is_target", scoredPsms);
			int massDiff = column(col, "mass_diff", scoredPsms);
			int logit = column(col, "class_logit", scoredPsms);
			int winner = column(col, "is_spectrum_winner", scoredPsms);
			int peak = column(col, "peak_hypothesis_accepted", scoredPsms);
			// Written by the closed-then-open identification only.
			int stage = col.getOrDefault("stage", -1);
			int closedAccepted = col.getOrDefault("closed_accepted", -1);
			int last = Math.max(Math.max(Math.max(file, peptide), Math.max(target, massDiff)),
					Math.max(Math.max(logit, winner), Math.max(peak, Math.max(stage, closedAccepted))));
			String line;
			while ((line = br.readLine()) != null) {
				String[] cs = line.split("\t", -1);
				if (cs.length <= last || !"true".equalsIgnoreCase(cs[winner])) {
					continue;
				}
				double score = parse(cs[logit]);
				if (!Double.isFinite(score)) {
					throw new IOException("Genome test: spectrum winner without class_logit in " + scoredPsms
							+ " (" + cs[file] + ", " + cs[peptide] + "); the test needs one score scale");
				}
				boolean isTarget = "true".equalsIgnoreCase(cs[target]);
				if (stage >= 0 && "closed".equals(cs[stage])) {
					hits.add(cs[file], cs[peptide], isTarget, GenomeEvidenceTest.PeptideClass.CLOSED, score);
				} else if (closedAccepted >= 0 && "true".equalsIgnoreCase(cs[closedAccepted])) {
					// The closed stage identified this spectrum; its open winner is no identification,
					// but an open decoy still belongs to the pi_g window.
					if (!isTarget) {
						hits.addWindowDecoy(cs[file], cs[peptide], score);
					}
				} else {
					hits.add(cs[file], cs[peptide], isTarget,
							classify(parse(cs[massDiff]), "true".equalsIgnoreCase(cs[peak])), score);
				}
			}
		}
		return hits;
	}

	static GenomeEvidenceTest.PeptideClass classify(double massDiff, boolean onAcceptedPeak) {
		if (Double.isFinite(massDiff)) {
			long k = Math.round(massDiff / ISOTOPE_SPACING);
			if (k >= 0 && k <= MAX_ISOTOPE_ERROR
					&& Math.abs(massDiff - k * ISOTOPE_SPACING) < UNMODIFIED_TOLERANCE_DA) {
				return GenomeEvidenceTest.PeptideClass.UNMODIFIED;
			}
		}
		return onAcceptedPeak ? GenomeEvidenceTest.PeptideClass.KNOWN_SHIFT
				: GenomeEvidenceTest.PeptideClass.UNEXPLAINED;
	}

	private static Set<String> genomesOf(File fasta, String decoyPrefix) throws IOException {
		Set<String> genomes = new HashSet<>();
		try (BufferedReader br = Files.newBufferedReader(fasta.toPath(), StandardCharsets.UTF_8)) {
			String line;
			while ((line = br.readLine()) != null) {
				if (line.startsWith(">")) {
					String accession = accession(line);
					if (!accession.startsWith(decoyPrefix)) {
						int idx = accession.indexOf('_');
						genomes.add(idx > 0 ? accession.substring(0, idx) : accession);
					}
				}
			}
		}
		return genomes;
	}

	/** Each accession once across the runs' databases, which share most of their genomes. */
	private static void forEachProtein(List<File> fastas, GenomeEvidenceTest.ProteinConsumer consumer)
			throws IOException {
		Set<String> seen = new HashSet<>();
		for (File fasta : fastas) {
			try (BufferedReader br = Files.newBufferedReader(fasta.toPath(), StandardCharsets.UTF_8)) {
				String name = null;
				StringBuilder seq = new StringBuilder();
				String line;
				while ((line = br.readLine()) != null) {
					if (line.startsWith(">")) {
						if (name != null && seq.length() > 0) {
							consumer.accept(name, seq.toString());
						}
						String accession = accession(line);
						name = seen.add(accession) ? accession : null;
						seq.setLength(0);
					} else if (name != null) {
						seq.append(line.trim());
					}
				}
				if (name != null && seq.length() > 0) {
					consumer.accept(name, seq.toString());
				}
			}
		}
	}

	private static String accession(String header) {
		String name = header.substring(1).trim();
		int sp = name.indexOf(' ');
		return sp > 0 ? name.substring(0, sp) : name;
	}

	private static int column(Map<String, Integer> col, String name, File file) throws IOException {
		Integer i = col.get(name);
		if (i == null) {
			throw new IOException("Genome test: column " + name + " missing from " + file);
		}
		return i;
	}

	private static double parse(String s) {
		if (s == null || s.isEmpty()) {
			return Double.NaN;
		}
		try {
			return Double.parseDouble(s);
		} catch (NumberFormatException e) {
			return Double.NaN;
		}
	}
}
