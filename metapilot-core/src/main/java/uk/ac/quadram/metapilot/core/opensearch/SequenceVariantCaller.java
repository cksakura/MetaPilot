package uk.ac.quadram.metapilot.core.opensearch;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.IdentityHashMap;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.Set;
import java.util.TreeSet;
import java.util.function.Consumer;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

/**
 * Calls reproducible single-residue sequence variants from open-search PSMs
 * whose delta mass {@link SubstitutionReassigner} has already explained as an
 * amino-acid substitution to a sequence absent from every searched database.
 *
 * <h2>What a call means</h2>
 * A variant here is a peptide that (a) is not in any database this run
 * searched, (b) differs from a database peptide by one residue at the position
 * MSFragger localized the delta to, and (c) is observed reproducibly. It is a
 * statement about the <em>sample's proteome relative to the nearest database
 * sequence</em>. It does not say how the difference arose: a strain that
 * differs from its reference genome and a mutation that arose in culture look
 * identical to the mass spectrometer. Only sequencing the sample separates them.
 *
 * <h2>Why these criteria</h2>
 * Measured on 12mix, run opensearch_MAG_search_0918/combined_msfragger_sage
 * (2026-09-28), on the 1,184 target CH2 (+/-14.0157 Da) PSMs whose substituted
 * candidate was absent from all 32 selected genomes, grouped by site and by how
 * many of the 9 files the same site was seen in:
 * <pre>
 *   files   sites  PSMs  reference peptide ever   two-nucleotide (N&lt;-&gt;Q)
 *                        seen unmodified
 *     1      319    321        54%                      12.9%
 *     2       88    183        43%                      10.2%
 *   3-5       87    342        37%                       3.4%
 *   6-9       37    338        22%                       2.7%
 * </pre>
 * One-off sites carry four times the two-nucleotide rate - a change no single
 * point mutation produces - and usually co-occur with the unmodified reference:
 * search artifacts or modifications. Recurrent sites are single-nucleotide,
 * reproducible, and usually replace the reference outright. Hence three pieces of
 * evidence, reported separately rather than folded into one score:
 * <ol>
 * <li><b>Recurrence</b> - distinct samples carrying the same site.</li>
 * <li><b>Codon reachability</b> - one nucleotide change under the standard
 * genetic code. Not a filter on real biology (multi-nucleotide differences exist
 * between strains) but the strongest single artifact signal measured.</li>
 * <li><b>Reference observation</b> - whether the database peptide is also seen
 * unmodified. Absent: the sample carries the variant instead of the reference.
 * Present: a mixed population, a reference peptide shared with another genome,
 * or a sub-stoichiometric modification - the three are not separable here.</li>
 * </ol>
 * A decoy control runs the same calls on decoy PSMs; on the run above it gave 2
 * decoy SEQUENCE_VARIANT calls against 272 target calls.
 *
 * <h2>What entrapment does and does not say about a call</h2>
 * On the same run 13.6% of SEQUENCE_VARIANT calls had a reference peptide found
 * only in non-spiked genomes. That is not an error rate: the reference is by
 * definition the nearest database sequence, and when the sample's protein is
 * absent from its own genome the nearest match is often a relative. The genome
 * of the reference is reported so this can be read per call.
 *
 * <h2>Isobaric modifications</h2>
 * Some substitutions have exactly the mass of a modification on the same
 * residue: D-&gt;E and methyl-D (methyl ester), G-&gt;A and methyl-G, S-&gt;T
 * and methyl-S, N-&gt;Q and methyl-N; N-&gt;D and deamidation. Backbone fragments
 * cannot separate them. They are listed per call in
 * {@code isobaric_modifications}; recurrence with the reference absent argues
 * for the variant (a sub-stoichiometric modification leaves the reference
 * visible), but a modification at full occupancy would look the same.
 */
public final class SequenceVariantCaller {

	private static final Logger LOGGER = LogManager.getLogger(SequenceVariantCaller.class);

	/** Evidence class of one variant. Ordered from strongest to weakest. */
	public enum VariantClass {
		/** Recurrent, single-nucleotide, reference peptide never observed unmodified. */
		SEQUENCE_VARIANT,
		/** Recurrent, single-nucleotide, reference peptide also observed unmodified. */
		SEQUENCE_VARIANT_WITH_REFERENCE,
		/** Single-nucleotide but seen in fewer than {@link Config#minSamples} samples. */
		LOW_RECURRENCE,
		/** Needs two or more nucleotide changes; the dominant artifact signature. */
		MULTI_NUCLEOTIDE
	}

	public static final class Config {
		public boolean enabled = true;
		/**
		 * Distinct samples a site must be seen in. Capped at the number of samples
		 * in the run, so a single-sample run still reports classes - but with no
		 * recurrence evidence behind them, which the table header states.
		 */
		public int minSamples = 3;
		/** |normalized delta| below this counts as an unmodified observation. */
		public double zeroWindowDa = 0.05;
		/** Tolerance for listing a modification as isobaric with the substitution. */
		public double isobaricToleranceDa = 0.005;
		/** Decoy accession prefix, skipped when reading protein sequences. */
		public String decoyTag = "rev_";
	}

	/** One called site, possibly supported by several overlapping peptides. */
	public static final class Variant {
		public String id;
		public VariantClass variantClass;
		public char from;
		public char to;
		public double substitutionMassDa;
		public boolean singleNucleotide;
		public final Set<String> isobaricModifications = new TreeSet<>();
		/** {@code accession:F123T}, 1-based protein position. */
		public final Set<String> proteinSites = new TreeSet<>();
		public final Set<String> referenceGenomes = new TreeSet<>();
		public final Set<String> referencePeptides = new TreeSet<>();
		public final Set<String> variantPeptides = new TreeSet<>();
		/** {@code REFERENCEPEPTIDE:F7T}, 1-based peptide position. */
		public final Set<String> peptideSites = new TreeSet<>();
		public final Set<String> samples = new TreeSet<>();
		public final Set<String> spectra = new HashSet<>();
		public final Set<Integer> peakIds = new TreeSet<>();
		public final List<Double> deltas = new ArrayList<>();
		public final List<OpenSearchPsm> psms = new ArrayList<>();
		public double bestQValue = Double.POSITIVE_INFINITY;
		/** Whether any supporting peptide carries the change at its N-terminal residue. */
		public boolean peptideNTerm;
		/** Samples in which any reference peptide of this site was seen unmodified. */
		public final Set<String> referenceSamples = new TreeSet<>();

		public int coObservedSamples() {
			int n = 0;
			for (String s : samples) {
				if (referenceSamples.contains(s)) {
					n++;
				}
			}
			return n;
		}

		public double medianDelta() {
			if (deltas.isEmpty()) {
				return Double.NaN;
			}
			List<Double> s = new ArrayList<>(deltas);
			Collections.sort(s);
			int m = s.size() / 2;
			return s.size() % 2 == 1 ? s.get(m) : 0.5 * (s.get(m - 1) + s.get(m));
		}

		/** Residue as a mass spectrometer reports it: I and L cannot be told apart. */
		public String substitution() {
			return residueLabel(from) + "->" + residueLabel(to);
		}
	}

	public static final class Result {
		public final List<Variant> variants;
		public final Map<OpenSearchPsm, Variant> byPsm;
		public final int samplesTotal;
		public final int effectiveMinSamples;
		public final int evidencePsms;
		public final Map<VariantClass, Integer> decoyCalls;
		public final boolean proteinPositionsAvailable;

		Result(List<Variant> variants, Map<OpenSearchPsm, Variant> byPsm, int samplesTotal, int effectiveMinSamples,
				int evidencePsms, Map<VariantClass, Integer> decoyCalls, boolean proteinPositionsAvailable) {
			this.variants = variants;
			this.byPsm = byPsm;
			this.samplesTotal = samplesTotal;
			this.effectiveMinSamples = effectiveMinSamples;
			this.evidencePsms = evidencePsms;
			this.decoyCalls = decoyCalls;
			this.proteinPositionsAvailable = proteinPositionsAvailable;
		}

		public int count(VariantClass c) {
			int n = 0;
			for (Variant v : variants) {
				if (v.variantClass == c) {
					n++;
				}
			}
			return n;
		}

		public static Result empty() {
			return new Result(new ArrayList<>(), new IdentityHashMap<>(), 0, 0, 0, new LinkedHashMap<>(), false);
		}
	}

	private SequenceVariantCaller() {
	}

	/**
	 * @param validated     the validated PSM set, targets and decoys; decoys drive
	 *                      the control count only
	 * @param proteinFastas databases this run searched, for protein positions;
	 *                      empty or missing files degrade to peptide-level sites
	 */
	public static Result call(List<OpenSearchPsm> validated, List<File> proteinFastas, Config cfg,
			Consumer<String> progress) {
		if (cfg == null) {
			cfg = new Config();
		}
		if (!cfg.enabled || validated == null || validated.isEmpty()) {
			return Result.empty();
		}

		Set<String> allSamples = new TreeSet<>();
		// Unmodified observations of each peptide, I/L-normalised, targets only.
		Map<String, Set<String>> unmodifiedSamples = new HashMap<>();
		List<OpenSearchPsm> targetEvidence = new ArrayList<>();
		List<OpenSearchPsm> decoyEvidence = new ArrayList<>();
		for (OpenSearchPsm psm : validated) {
			allSamples.add(psm.getFileName());
			if (psm.isTarget() && Math.abs(psm.getNormalizedMassDiff()) < cfg.zeroWindowDa
					&& psm.getPeptide() != null) {
				unmodifiedSamples.computeIfAbsent(normalise(psm.getPeptide()), k -> new HashSet<>())
						.add(psm.getFileName());
			}
			if (!SubstitutionReassigner.Verdict.SUBSTITUTION_NOT_IN_DB.name().equals(psm.getSubstitutionVerdict())
					|| siteOf(psm) == null) {
				continue;
			}
			(psm.isTarget() ? targetEvidence : decoyEvidence).add(psm);
		}
		int minSamples = Math.max(1, Math.min(cfg.minSamples, allSamples.size()));

		// Protein sequences for every accession a target evidence PSM names.
		Set<String> wanted = new HashSet<>();
		for (OpenSearchPsm psm : targetEvidence) {
			wanted.addAll(proteinsOf(psm));
		}
		Map<String, String> sequences = loadSequences(proteinFastas, wanted, cfg.decoyTag);

		List<Variant> variants = group(targetEvidence, sequences);
		for (Variant v : variants) {
			for (String ref : v.referencePeptides) {
				Set<String> seen = unmodifiedSamples.get(normalise(ref));
				if (seen != null) {
					v.referenceSamples.addAll(seen);
				}
			}
			classify(v, minSamples, cfg);
		}
		variants.sort((a, b) -> {
			int c = a.variantClass.compareTo(b.variantClass);
			if (c != 0) {
				return c;
			}
			c = Integer.compare(b.samples.size(), a.samples.size());
			return c != 0 ? c : Integer.compare(b.spectra.size(), a.spectra.size());
		});
		Map<OpenSearchPsm, Variant> byPsm = new IdentityHashMap<>();
		for (int i = 0; i < variants.size(); i++) {
			Variant v = variants.get(i);
			v.id = String.format(Locale.ROOT, "SV%05d", i + 1);
			for (OpenSearchPsm psm : v.psms) {
				byPsm.put(psm, v);
			}
		}

		// Decoy control: the same grouping and classes on decoy PSMs. No protein
		// positions and no reference observations - a decoy's reference is a
		// decoy and is never observed unmodified as a target.
		Map<VariantClass, Integer> decoyCalls = new LinkedHashMap<>();
		for (VariantClass c : VariantClass.values()) {
			decoyCalls.put(c, 0);
		}
		for (Variant v : group(decoyEvidence, Collections.emptyMap())) {
			classify(v, minSamples, cfg);
			decoyCalls.merge(v.variantClass, 1, Integer::sum);
		}

		Result result = new Result(variants, byPsm, allSamples.size(), minSamples, targetEvidence.size(), decoyCalls,
				!sequences.isEmpty());
		String line = String.format(Locale.ROOT,
				"Sequence variants: %,d substitution PSMs absent from the searched databases -> %,d sites; "
						+ "%d SEQUENCE_VARIANT, %d SEQUENCE_VARIANT_WITH_REFERENCE, %d LOW_RECURRENCE, "
						+ "%d MULTI_NUCLEOTIDE (recurrence >= %d of %d samples). Decoy control: %d SEQUENCE_VARIANT, "
						+ "%d SEQUENCE_VARIANT_WITH_REFERENCE",
				targetEvidence.size(), variants.size(), result.count(VariantClass.SEQUENCE_VARIANT),
				result.count(VariantClass.SEQUENCE_VARIANT_WITH_REFERENCE),
				result.count(VariantClass.LOW_RECURRENCE), result.count(VariantClass.MULTI_NUCLEOTIDE), minSamples,
				allSamples.size(), decoyCalls.get(VariantClass.SEQUENCE_VARIANT),
				decoyCalls.get(VariantClass.SEQUENCE_VARIANT_WITH_REFERENCE));
		if (progress != null) {
			progress.accept(line);
			if (sequences.isEmpty() && !targetEvidence.isEmpty()) {
				progress.accept("Sequence variants: no protein sequences could be read, so sites are reported "
						+ "per peptide and overlapping peptides of one site are not merged.");
			}
		}
		LOGGER.info(line);
		return result;
	}

	private static void classify(Variant v, int minSamples, Config cfg) {
		v.singleNucleotide = GeneticCode.singleNucleotide(v.from, v.to);
		v.substitutionMassDa = SubstitutionMassTable.residueMass(v.to) - SubstitutionMassTable.residueMass(v.from);
		v.isobaricModifications.addAll(IsobaricModificationTable.sameResidue(v.from, v.substitutionMassDa,
				v.peptideNTerm, cfg.isobaricToleranceDa));
		if (!v.singleNucleotide) {
			v.variantClass = VariantClass.MULTI_NUCLEOTIDE;
		} else if (v.samples.size() < minSamples) {
			v.variantClass = VariantClass.LOW_RECURRENCE;
		} else if (v.referenceSamples.isEmpty()) {
			v.variantClass = VariantClass.SEQUENCE_VARIANT;
		} else {
			v.variantClass = VariantClass.SEQUENCE_VARIANT_WITH_REFERENCE;
		}
	}

	/**
	 * Groups evidence PSMs into sites. With protein sequences, two peptides that
	 * put the same residue change on the same protein position are one site -
	 * this is what merges missed-cleavage forms. Without, the site is the
	 * (reference, variant) peptide pair.
	 */
	private static List<Variant> group(List<OpenSearchPsm> evidence, Map<String, String> sequences) {
		// Union-find over evidence PSMs, joined through shared site tokens.
		int n = evidence.size();
		int[] parent = new int[n];
		for (int i = 0; i < n; i++) {
			parent[i] = i;
		}
		Map<String, Integer> firstByToken = new HashMap<>();
		List<List<String>> proteinSitesOf = new ArrayList<>(n);
		for (int i = 0; i < n; i++) {
			OpenSearchPsm psm = evidence.get(i);
			Site site = siteOf(psm);
			List<String> tokens = new ArrayList<>();
			List<String> proteinSites = new ArrayList<>();
			for (String acc : proteinsOf(psm)) {
				String seq = sequences.get(acc);
				if (seq == null) {
					continue;
				}
				for (int start : occurrences(seq, psm.getPeptide())) {
					int pos = start + site.position + 1;
					String token = acc + ":" + residueLabel(site.from) + pos + residueLabel(site.to);
					proteinSites.add(token);
					tokens.add(token);
				}
			}
			if (tokens.isEmpty()) {
				tokens.add("pep:" + normalise(psm.getPeptide()) + ">" + normalise(site.candidate));
			}
			proteinSitesOf.add(proteinSites);
			for (String t : tokens) {
				Integer j = firstByToken.putIfAbsent(t, i);
				if (j != null) {
					union(parent, i, j);
				}
			}
		}
		Map<Integer, Variant> byRoot = new LinkedHashMap<>();
		for (int i = 0; i < n; i++) {
			OpenSearchPsm psm = evidence.get(i);
			Site site = siteOf(psm);
			Variant v = byRoot.computeIfAbsent(find(parent, i), k -> new Variant());
			v.from = site.from;
			v.to = site.to;
			v.peptideNTerm |= site.position == 0;
			v.psms.add(psm);
			v.referencePeptides.add(psm.getPeptide());
			v.variantPeptides.add(site.candidate);
			v.peptideSites.add(psm.getPeptide() + ":" + residueLabel(site.from) + (site.position + 1)
					+ residueLabel(site.to));
			v.proteinSites.addAll(proteinSitesOf.get(i));
			for (String acc : proteinsOf(psm)) {
				v.referenceGenomes.add(genomeOf(acc));
			}
			v.samples.add(psm.getFileName());
			v.spectra.add(psm.getFileName() + "#" + psm.getScan());
			if (psm.getPeakId() >= 0) {
				v.peakIds.add(psm.getPeakId());
			}
			v.deltas.add(psm.getNormalizedMassDiff());
			v.bestQValue = Math.min(v.bestQValue, psm.getQValue());
		}
		return new ArrayList<>(byRoot.values());
	}

	private static int find(int[] parent, int i) {
		while (parent[i] != i) {
			parent[i] = parent[parent[i]];
			i = parent[i];
		}
		return i;
	}

	private static void union(int[] parent, int a, int b) {
		int ra = find(parent, a);
		int rb = find(parent, b);
		if (ra != rb) {
			parent[Math.max(ra, rb)] = Math.min(ra, rb);
		}
	}

	/** The single residue change between a PSM's peptide and its candidate. */
	private static final class Site {
		final int position;
		final char from;
		final char to;
		final String candidate;

		Site(int position, char from, char to, String candidate) {
			this.position = position;
			this.from = from;
			this.to = to;
			this.candidate = candidate;
		}
	}

	/**
	 * Derived by comparing the two sequences, not by parsing
	 * {@code deltaMassReassignment}'s text: the sequences are the fact, the label
	 * is a rendering of it. Null unless they differ at exactly one position.
	 */
	private static Site siteOf(OpenSearchPsm psm) {
		String p = psm.getPeptide();
		String c = psm.getSubstitutionCandidate();
		if (p == null || c == null || p.length() != c.length()) {
			return null;
		}
		int at = -1;
		for (int i = 0; i < p.length(); i++) {
			if (p.charAt(i) != c.charAt(i)) {
				if (at >= 0) {
					return null;
				}
				at = i;
			}
		}
		return at < 0 ? null : new Site(at, p.charAt(at), c.charAt(at), c);
	}

	private static List<String> proteinsOf(OpenSearchPsm psm) {
		List<String> out = new ArrayList<>();
		String[] all = psm.getAllProteins();
		if (all != null) {
			for (String a : all) {
				if (a != null && !a.isBlank()) {
					out.add(a.trim());
				}
			}
		}
		if (out.isEmpty() && psm.getProtein() != null && !psm.getProtein().isBlank()) {
			out.add(psm.getProtein().trim());
		}
		return out;
	}

	/** Start offsets of {@code peptide} in {@code protein}, I and L treated as equal. */
	private static List<Integer> occurrences(String protein, String peptide) {
		List<Integer> out = new ArrayList<>();
		String hay = normalise(protein);
		String needle = normalise(peptide);
		int from = 0;
		int at;
		while ((at = hay.indexOf(needle, from)) >= 0) {
			out.add(at);
			from = at + 1;
		}
		return out;
	}

	/** Streams the FASTAs once, keeping only the accessions asked for. */
	private static Map<String, String> loadSequences(List<File> fastas, Set<String> wanted, String decoyTag) {
		Map<String, String> out = new HashMap<>();
		if (fastas == null || wanted.isEmpty()) {
			return out;
		}
		String tag = decoyTag == null ? "" : decoyTag.toLowerCase(Locale.ROOT);
		for (File fasta : fastas) {
			if (fasta == null || !fasta.isFile() || out.size() == wanted.size()) {
				continue;
			}
			try (BufferedReader br = new BufferedReader(new FileReader(fasta))) {
				String line;
				String name = null;
				StringBuilder seq = new StringBuilder();
				while ((line = br.readLine()) != null) {
					if (line.startsWith(">")) {
						if (name != null) {
							out.put(name, seq.toString());
						}
						seq.setLength(0);
						String h = line.substring(1).trim();
						int sp = h.indexOf(' ');
						h = sp > 0 ? h.substring(0, sp) : h;
						boolean keep = wanted.contains(h) && !out.containsKey(h)
								&& (tag.isEmpty() || !h.toLowerCase(Locale.ROOT).startsWith(tag));
						name = keep ? h : null;
					} else if (name != null) {
						seq.append(line.trim());
					}
				}
				if (name != null) {
					out.put(name, seq.toString());
				}
			} catch (IOException e) {
				LOGGER.warn("Sequence variants: could not read " + fasta + "; protein positions from it are missing",
						e);
			}
		}
		return out;
	}

	private static final Pattern MAG_ACCESSION = Pattern.compile("^(.+)_\\d+$");

	/**
	 * Genome of a MAG protein accession ({@code MGYG000002463_01761} ->
	 * {@code MGYG000002463}). Any other accession is returned unchanged.
	 */
	static String genomeOf(String accession) {
		Matcher m = MAG_ACCESSION.matcher(accession);
		return m.matches() ? m.group(1) : accession;
	}

	private static String normalise(String s) {
		return s.replace('I', 'L');
	}

	private static String residueLabel(char aa) {
		return aa == 'I' || aa == 'L' ? "I/L" : String.valueOf(aa);
	}

	/** The change part of a protein-site token, e.g. {@code T123S} or {@code V45I/L}. */
	private static final Pattern PROTEIN_SITE = Pattern.compile("^(I/L|[A-Z])(\\d+)(I/L|[A-Z])$");

	/** See {@link IsobaricModificationTable#appliesTo}. */
	static boolean appliesTo(String title, char residue, boolean peptideNTerm) {
		return IsobaricModificationTable.appliesTo(title, residue, peptideNTerm);
	}

	/** Standard genetic code; used only to ask whether one nucleotide change suffices. */
	static final class GeneticCode {
		private static final String BASES = "TCAG";
		private static final String AA = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG";
		private static final Map<String, Character> CODE = new HashMap<>();
		static {
			int k = 0;
			for (char a : BASES.toCharArray()) {
				for (char b : BASES.toCharArray()) {
					for (char c : BASES.toCharArray()) {
						CODE.put("" + a + b + c, AA.charAt(k++));
					}
				}
			}
		}

		private GeneticCode() {
		}

		/**
		 * True when some codon of {@code from} becomes a codon of {@code to} by one
		 * substitution. I and L count as one residue: a spectrum cannot tell them
		 * apart, so either reading that is reachable is enough.
		 */
		static boolean singleNucleotide(char from, char to) {
			for (char f : isoforms(from)) {
				for (char t : isoforms(to)) {
					if (reachable(f, t)) {
						return true;
					}
				}
			}
			return false;
		}

		private static char[] isoforms(char aa) {
			return aa == 'I' || aa == 'L' ? new char[] { 'I', 'L' } : new char[] { aa };
		}

		private static boolean reachable(char from, char to) {
			for (Map.Entry<String, Character> e : CODE.entrySet()) {
				if (e.getValue() != from) {
					continue;
				}
				char[] codon = e.getKey().toCharArray();
				for (int i = 0; i < 3; i++) {
					char orig = codon[i];
					for (char b : BASES.toCharArray()) {
						if (b == orig) {
							continue;
						}
						codon[i] = b;
						Character aa = CODE.get(new String(codon));
						if (aa != null && aa == to) {
							return true;
						}
					}
					codon[i] = orig;
				}
			}
			return false;
		}
	}

	// ---------------------------------------------------------------- output

	/** Site-level table: one row per called variant. */
	public static void writeVariantTable(Result result, Config cfg, File out) throws IOException {
		out.getParentFile().mkdirs();
		try (PrintWriter w = new PrintWriter(out, "UTF-8")) {
			w.println("# Single-residue sequence variants: peptides absent from every searched database that differ");
			w.println("# from a database peptide by one residue at the position the open-search delta localized to.");
			w.println("# A call describes the sample relative to the nearest database sequence; it cannot say");
			w.println("# whether the difference is a strain-vs-reference difference or a mutation that arose in");
			w.println("# culture. I and L are indistinguishable by mass and are written I/L.");
			w.println("#");
			w.println("# variant_class:");
			w.println("#   SEQUENCE_VARIANT                 seen in >= min_samples samples, one-nucleotide change,");
			w.println("#                                    reference peptide never observed unmodified");
			w.println("#   SEQUENCE_VARIANT_WITH_REFERENCE  as above, but the reference peptide is also observed:");
			w.println("#                                    mixed population, reference shared with another genome,");
			w.println("#                                    or a sub-stoichiometric modification");
			w.println("#   LOW_RECURRENCE                   one-nucleotide change seen in < min_samples samples");
			w.println("#   MULTI_NUCLEOTIDE                 needs >= 2 nucleotide changes (dominant artifact signature)");
			w.println("#");
			w.println("# min_samples=" + result.effectiveMinSamples + " (configured " + cfg.minSamples + ")"
					+ "  samples_in_run=" + result.samplesTotal + "  evidence_psms=" + result.evidencePsms
					+ "  protein_positions=" + (result.proteinPositionsAvailable ? "yes" : "no"));
			if (result.samplesTotal < cfg.minSamples) {
				w.println("# WARNING: fewer samples than min_samples; recurrence is not evidence in this run.");
			}
			StringBuilder decoy = new StringBuilder("# decoy_control:");
			for (Map.Entry<VariantClass, Integer> e : result.decoyCalls.entrySet()) {
				decoy.append(' ').append(e.getKey().name()).append('=').append(e.getValue());
			}
			w.println(decoy);
			w.println("variant_id\tvariant_class\tsubstitution\tsubstitution_mass_da\tmedian_delta_mass_da\t"
					+ "single_nucleotide\tisobaric_modifications\tprotein_sites\treference_genomes\t"
					+ "reference_peptides\tvariant_peptides\tpeptide_sites\tsamples_observed\tsamples_total\t"
					+ "spectra\tpsms\tbest_q_value\treference_observed_samples\treference_co_observed_samples\t"
					+ "samples\tpeak_ids");
			for (Variant v : result.variants) {
				List<String> peaks = new ArrayList<>();
				for (Integer p : v.peakIds) {
					peaks.add(String.valueOf(p));
				}
				w.println(v.id + "\t" + v.variantClass + "\t" + v.substitution() + "\t"
						+ String.format(Locale.ROOT, "%.5f", v.substitutionMassDa) + "\t"
						+ String.format(Locale.ROOT, "%.5f", v.medianDelta()) + "\t" + v.singleNucleotide + "\t"
						+ join(v.isobaricModifications) + "\t" + join(v.proteinSites) + "\t"
						+ join(v.referenceGenomes) + "\t" + join(v.referencePeptides) + "\t"
						+ join(v.variantPeptides) + "\t" + join(v.peptideSites) + "\t" + v.samples.size() + "\t"
						+ result.samplesTotal + "\t" + v.spectra.size() + "\t" + v.psms.size() + "\t"
						+ String.format(Locale.ROOT, "%.6g", v.bestQValue) + "\t" + v.referenceSamples.size() + "\t"
						+ v.coObservedSamples() + "\t" + join(v.samples) + "\t" + String.join(";", peaks));
			}
		}
	}

	/**
	 * Variant proteins for the recurrent classes: each reference protein with the
	 * one residue replaced, for a follow-up closed search or quantification that
	 * should see the variant peptide as an ordinary unmodified peptide. For an
	 * I/L site the residue written is the one the variant peptide carries.
	 *
	 * @return number of protein entries written
	 */
	public static int writeVariantFasta(Result result, List<File> proteinFastas, Config cfg, File out)
			throws IOException {
		Set<String> wanted = new HashSet<>();
		List<Variant> recurrent = new ArrayList<>();
		for (Variant v : result.variants) {
			if (v.variantClass == VariantClass.SEQUENCE_VARIANT
					|| v.variantClass == VariantClass.SEQUENCE_VARIANT_WITH_REFERENCE) {
				recurrent.add(v);
				for (String site : v.proteinSites) {
					wanted.add(site.substring(0, site.lastIndexOf(':')));
				}
			}
		}
		Map<String, String> sequences = loadSequences(proteinFastas, wanted, cfg.decoyTag);
		int written = 0;
		out.getParentFile().mkdirs();
		try (PrintWriter w = new PrintWriter(out, "UTF-8")) {
			for (Variant v : recurrent) {
				for (String site : v.proteinSites) {
					int colon = site.lastIndexOf(':');
					String acc = site.substring(0, colon);
					String seq = sequences.get(acc);
					Matcher m = PROTEIN_SITE.matcher(site.substring(colon + 1));
					if (seq == null || !m.matches()) {
						continue;
					}
					int pos = Integer.parseInt(m.group(2)) - 1;
					if (pos < 0 || pos >= seq.length()) {
						continue;
					}
					String variantSeq = seq.substring(0, pos) + v.to + seq.substring(pos + 1);
					w.println(">" + acc + "|" + v.id + "|" + site.substring(colon + 1) + " class=" + v.variantClass
							+ " samples=" + v.samples.size() + "/" + result.samplesTotal);
					for (int i = 0; i < variantSeq.length(); i += 60) {
						w.println(variantSeq.substring(i, Math.min(variantSeq.length(), i + 60)));
					}
					written++;
				}
			}
		}
		return written;
	}

	private static String join(Set<String> s) {
		return String.join(";", s);
	}
}
