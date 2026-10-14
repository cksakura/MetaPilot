package uk.ac.quadram.metapilot.core.opensearch;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.Set;
import java.util.function.Consumer;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

/**
 * Is a substitution candidate's partner genome actually plausible for the
 * sample it was found in - not just "somewhere in the searched database"?
 *
 * <h2>The problem this exists for</h2>
 * {@link SubstitutionReassigner} and {@link TrypticPeptideIndex} answer "does
 * this substituted sequence exist in the searched database", which is genome-
 * blind (see {@code TrypticPeptideIndex}'s own javadoc). Resolved to the
 * genome level on the 12mix run (2026-09-16/17), a database match on its own
 * carries a real error rate: of 899 matched original-to-form peptide pairs
 * whose form landed in a genome the original was not already in, only 63.6%
 * had that partner be one of the twelve spiked (true) genomes - 36.4% were
 * entrapment. Splitting the same pairs by the best-supported partner genome's
 * own {@code TotalUniqueMarkers} (the same statistic stage-1 genome selection
 * already computes, in {@code <sample>/marker/marker_genome_ranked.tsv})
 * separates the two almost perfectly:
 * <pre>
 *   max unique markers of partner genome   true   entrapment   true %
 *   0                                         0          108     0.0%
 *   1-4                                     101          194    34.2%
 *   5-19                                    325           23    93.4%
 *   &gt;=20                                    142            0   100.0%
 * </pre>
 * As a threshold on the same data:
 * <pre>
 *   threshold   TP    FP    precision   recall
 *   &gt;=0 (none)  568   325    63.6%      100%
 *   &gt;=1         568   217    72.4%      100%
 *   &gt;=3         511    89    85.2%       90.0%
 *   &gt;=5         467    23    95.3%       82.2%
 *   &gt;=8         379    12    96.9%       66.7%
 *   &gt;=20        142     0   100.0%       25.0%
 * </pre>
 * {@code TotalPSM} does NOT separate them (median 2,194 true vs 1,744
 * entrapment on the same pairs) - entrapment genomes are close relatives and
 * inherit shared peptides, so only unique-marker evidence is discriminative.
 *
 * <h2>What this deliberately does not answer</h2>
 * <ul>
 * <li>This tests "is the partner genome present in this sample", not "is this
 * sequence a plausible variant of the assigned protein". They are different
 * questions; an unrelated protein in a well-supported genome still passes
 * here.</li>
 * <li>It inherits any bias in stage-1 genome selection. A true genome that
 * selection missed cannot supply a plausible partner - this is a filter on
 * partners, never on the underlying spectrum.</li>
 * </ul>
 * A complementary, weaker statistic - how many genomes carry the candidate
 * form at all ("partner conservation") - is deliberately not folded in here:
 * a form present in most of the searched database is uninformative about
 * strain identity on its own, so it is reported alongside this class's
 * verdict rather than combined with it.
 *
 * <h2>Per sample, not pooled</h2>
 * The threshold table above was computed on the median across samples; this
 * class answers per sample instead, which is strictly better - a genome can
 * be well supported in one sample and absent from another, and pooling would
 * hide that.
 */
public final class PartnerGenomeSupport {

	private static final Logger LOGGER = LogManager.getLogger(PartnerGenomeSupport.class);

	/** Filename stage-1 genome selection writes under each sample's marker directory. */
	public static final String MARKER_GENOME_RANKED_FILE = "marker_genome_ranked.tsv";

	public static final class Config {
		/**
		 * A partner genome is plausible when its {@code TotalUniqueMarkers} in the
		 * sample meets or exceeds this. Measured on 12mix (see class javadoc):
		 * 95.3% precision at 82.2% recall against the twelve spiked genomes. 3 is
		 * the reasonable recall-favouring alternative (85.2% / 90.0%); 0 disables
		 * the test and restores "present in the database" as sufficient, with its
		 * measured 36.4% error rate.
		 */
		public int minPartnerUniqueMarkers = 5;
	}

	private final Map<String, Map<String, Integer>> uniqueMarkersBySampleAndGenome;
	private final int samplesLoaded;

	private PartnerGenomeSupport(Map<String, Map<String, Integer>> data, int samplesLoaded) {
		this.uniqueMarkersBySampleAndGenome = data;
		this.samplesLoaded = samplesLoaded;
	}

	public int samplesLoaded() {
		return samplesLoaded;
	}

	/**
	 * {@code TotalUniqueMarkers} for this genome in this sample, or 0 when the
	 * genome was never ranked there - stage-1 genome selection never considered
	 * it a candidate for this sample at all, which is the strongest possible
	 * "not plausible" signal, not a missing-data gap.
	 */
	public int uniqueMarkers(String sample, String genome) {
		if (sample == null || genome == null) {
			return 0;
		}
		Map<String, Integer> bySample = uniqueMarkersBySampleAndGenome.get(sample);
		if (bySample == null) {
			return 0;
		}
		Integer v = bySample.get(genome);
		return v == null ? 0 : v;
	}

	/**
	 * Best (max) unique-marker support among the given candidate genomes, for one
	 * sample. This is "plausibility" in the sense the class javadoc measures it -
	 * a raw marker count, not a normalised probability - so callers compare it
	 * against {@link Config#minPartnerUniqueMarkers} themselves, or call
	 * {@link #isPlausible}.
	 */
	public double plausibility(String sample, Set<String> genomes) {
		if (genomes == null || genomes.isEmpty()) {
			return 0.0;
		}
		int best = 0;
		for (String genome : genomes) {
			int markers = uniqueMarkers(sample, genome);
			if (markers > best) {
				best = markers;
			}
		}
		return best;
	}

	/** True when the best-supported genome in the set clears the configured threshold. */
	public boolean isPlausible(String sample, Set<String> genomes, Config cfg) {
		int threshold = (cfg == null ? new Config() : cfg).minPartnerUniqueMarkers;
		return plausibility(sample, genomes) >= threshold;
	}

	/**
	 * Loads every sample's {@code marker_genome_ranked.tsv}.
	 * <p>
	 * A sample whose file is missing is not an error - it is recorded as zero
	 * genomes ranked, so every partner genome tested against it reports as
	 * unranked (0 markers, not plausible) rather than silently skipping the
	 * sample or throwing.
	 */
	public static PartnerGenomeSupport load(File resultRoot, List<String[]> experiments, Consumer<String> progress)
			throws IOException {
		Map<String, Map<String, Integer>> data = new HashMap<>();
		int filesRead = 0;
		int missing = 0;
		for (String[] experiment : experiments) {
			String sample = baseName(experiment[0]);
			File rankedTsv = new File(OpenSearchLayout.markerDir(resultRoot, sample), MARKER_GENOME_RANKED_FILE);
			Map<String, Integer> byGenome = new HashMap<>();
			if (rankedTsv.isFile()) {
				readRanked(rankedTsv, byGenome);
				filesRead++;
			} else {
				missing++;
			}
			data.put(sample, byGenome);
		}
		if (progress != null) {
			progress.accept(String.format(Locale.ROOT,
					"Partner-genome support: loaded marker rankings for %d of %d sample(s)%s", filesRead,
					experiments.size(), missing > 0
							? " (" + missing + " missing " + MARKER_GENOME_RANKED_FILE
									+ " - those samples report every partner as unranked)"
							: ""));
		}
		LOGGER.info("PartnerGenomeSupport: {} of {} sample marker ranking(s) loaded", filesRead, experiments.size());
		return new PartnerGenomeSupport(data, filesRead);
	}

	private static void readRanked(File tsv, Map<String, Integer> out) throws IOException {
		try (BufferedReader br = new BufferedReader(new FileReader(tsv))) {
			String header = br.readLine();
			if (header == null) {
				return;
			}
			String[] cols = header.split("\t", -1);
			int genomeIdx = -1;
			int markersIdx = -1;
			for (int i = 0; i < cols.length; i++) {
				if ("Genome".equals(cols[i])) {
					genomeIdx = i;
				} else if ("TotalUniqueMarkers".equals(cols[i])) {
					markersIdx = i;
				}
			}
			if (genomeIdx < 0 || markersIdx < 0) {
				throw new IOException(
						tsv + " has no Genome/TotalUniqueMarkers column - stage-1 output format changed?");
			}
			String line;
			while ((line = br.readLine()) != null) {
				if (line.isEmpty()) {
					continue;
				}
				String[] f = line.split("\t", -1);
				if (f.length <= Math.max(genomeIdx, markersIdx)) {
					continue;
				}
				try {
					// A genome should appear once per sample's ranking; max is a defensive
					// tie-break, not an expected merge.
					out.merge(f[genomeIdx], Integer.parseInt(f[markersIdx].trim()), Math::max);
				} catch (NumberFormatException e) {
					// Malformed row; skip rather than fail the whole load.
				}
			}
		}
	}

	private static String baseName(String path) {
		String name = new File(path).getName();
		int dot = name.lastIndexOf('.');
		return dot > 0 ? name.substring(0, dot) : name;
	}

	/**
	 * One (origin genome, partner genome) row of {@code genome_crosstalk.tsv} -
	 * see {@link #writeCrosstalkEvidence}. Built by the caller from
	 * {@link SubstitutionReassigner} proposals or an equivalent genome-resolved
	 * substitution export; this class supplies only the marker lookup and the
	 * plausibility verdict, per {@code opensearch_postvalidation_round3_plan.md}
	 * section 5's "report-only, change nothing" scope for this component.
	 */
	public static final class CrosstalkEvidence {
		public final String originGenome;
		public final String partnerGenome;
		public final int nPeptides;
		public final long nPsms;
		public final List<String> topSubstitutions;
		/** Best support for partnerGenome across every sample it was checked in. */
		public final int maxUniqueMarkers;
		public final int samplesChecked;
		public final int samplesPlausible;
		public final boolean partnerPlausible;

		public CrosstalkEvidence(String originGenome, String partnerGenome, int nPeptides, long nPsms,
				List<String> topSubstitutions, int maxUniqueMarkers, int samplesChecked, int samplesPlausible,
				boolean partnerPlausible) {
			this.originGenome = originGenome;
			this.partnerGenome = partnerGenome;
			this.nPeptides = nPeptides;
			this.nPsms = nPsms;
			this.topSubstitutions = topSubstitutions;
			this.maxUniqueMarkers = maxUniqueMarkers;
			this.samplesChecked = samplesChecked;
			this.samplesPlausible = samplesPlausible;
			this.partnerPlausible = partnerPlausible;
		}
	}

	/**
	 * Evaluates one (origin, partner) pair against every loaded sample and
	 * returns its evidence row. {@code samples} is normally every sample this
	 * pair's supporting PSMs came from; when that breakdown is not available
	 * (a pooled, not-per-PSM export) the caller may pass every loaded sample,
	 * which checks partner plausibility everywhere rather than only where the
	 * evidence actually occurred - conservative in one direction (a genome
	 * plausible somewhere it was not actually observed still counts as
	 * checked), and this is recorded in {@code samplesChecked} so the
	 * approximation is visible in the output rather than silent.
	 */
	public CrosstalkEvidence evaluate(String originGenome, String partnerGenome, int nPeptides, long nPsms,
			List<String> topSubstitutions, Set<String> samples, Config cfg) {
		int threshold = (cfg == null ? new Config() : cfg).minPartnerUniqueMarkers;
		int max = 0;
		int plausibleCount = 0;
		int checked = 0;
		for (String sample : samples) {
			checked++;
			int markers = uniqueMarkers(sample, partnerGenome);
			if (markers > max) {
				max = markers;
			}
			if (markers >= threshold) {
				plausibleCount++;
			}
		}
		return new CrosstalkEvidence(originGenome, partnerGenome, nPeptides, nPsms, topSubstitutions, max, checked,
				plausibleCount, max >= threshold);
	}

	/**
	 * Writes {@code genome_crosstalk.tsv} - which genomes are drawing peptides
	 * from which, per {@code opensearch_postvalidation_round3_plan.md} section 9.
	 */
	public static void writeCrosstalkEvidence(List<CrosstalkEvidence> rows, Config cfg, File out) throws IOException {
		int threshold = (cfg == null ? new Config() : cfg).minPartnerUniqueMarkers;
		out.getParentFile().mkdirs();
		try (PrintWriter w = new PrintWriter(out, "UTF-8")) {
			w.println("# One row per (origin genome, partner genome) pair a substitution candidate proposes moving");
			w.println("# a peptide between. partner_plausible is true when the partner genome's own");
			w.println("# TotalUniqueMarkers (stage-1 genome selection, per sample) reaches " + threshold
					+ " in at least one");
			w.println("# checked sample - measured on 12mix at 95.3% precision / 82.2% recall against the spiked");
			w.println("# genomes. This tests presence, not homology: a plausible partner can still be an unrelated");
			w.println("# protein in a real genome. See PartnerGenomeSupport's class javadoc for the full threshold");
			w.println("# table and known limits.");
			w.println("origin_genome\tpartner_genome\tn_peptides\tn_psms\ttop_substitutions\t"
					+ "max_partner_unique_markers\tsamples_checked\tsamples_plausible\tpartner_plausible");
			for (CrosstalkEvidence r : rows) {
				w.println(r.originGenome + "\t" + r.partnerGenome + "\t" + r.nPeptides + "\t" + r.nPsms + "\t"
						+ String.join(",", r.topSubstitutions) + "\t" + r.maxUniqueMarkers + "\t" + r.samplesChecked
						+ "\t" + r.samplesPlausible + "\t" + r.partnerPlausible);
			}
		}
	}
}
