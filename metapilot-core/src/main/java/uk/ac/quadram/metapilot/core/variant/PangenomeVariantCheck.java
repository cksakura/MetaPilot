package uk.ac.quadram.metapilot.core.variant;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.io.Writer;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.StandardCopyOption;
import java.util.ArrayList;
import java.util.BitSet;
import java.util.Collection;
import java.util.Comparator;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.function.Consumer;

import uk.ac.quadram.metapilot.core.opensearch.IsobaricModificationTable;

/**
 * Pangenome variant check: is each substitution an open search called a
 * sequence that exists in nature (UHGP-100), and more often than a matched
 * random substitution would be?
 * <p>
 * Report-only add-on. It reads what a workflow has already decided and writes
 * five new files; nothing it computes feeds back into FDR, peak acceptance or
 * the {@code SubstitutionReassigner} gate, and no existing output changes.
 * Both engines that can run it (pFind open search and MetaPilot open search)
 * hand it {@link VariantCandidate}s from their own front end and share
 * everything from here on, so lookup, controls, classes and file formats are
 * one code path.
 * <p>
 * A failure is logged and reported in the {@link Outcome}, never thrown: the
 * workflow's own results are complete before this runs.
 */
public final class PangenomeVariantCheck {

	public static final String CANDIDATES_FILE = "variant_candidates.tsv";
	public static final String SITES_FILE = "variant_sites.tsv";
	public static final String SUMMARY_FILE = "variant_pangenome_summary.tsv";
	public static final String SPECIES_FILE = "variant_species_summary.tsv";
	/** What every column of the four tables above means. */
	public static final String COLUMNS_FILE = "variant_columns.tsv";

	private static final String LOG_PREFIX = "Pangenome variant check: ";
	private static final String NA = "NA";

	/** Everything one run needs; the front end fills the first block. */
	public static final class Request {
		public VariantCandidate.Engine engine;
		public List<VariantCandidate> candidates = new ArrayList<>();
		/** Every identified base peptide of this run and engine; the control pool. */
		public Collection<String> identifiedBasePeptides = new ArrayList<>();
		/**
		 * Every protein database this run searched, across all samples; target
		 * entries only are used. {@code in_searched_db} is true when the variant
		 * is in any of them, not necessarily the one its own spectrum was
		 * searched against.
		 */
		public List<File> searchedFastas = new ArrayList<>();
		public String decoyPrefix = "";

		public PangenomeVariantSettings settings = new PangenomeVariantSettings();
		/** {@code MagDbItem.getCurrentFile()}; locates the default resource and genome metadata. */
		public File magVersionDir;
		/** Holds {@code genome_final_selected.tsv}, the DETECTED species scope. */
		public File taxonomyDir;
		/**
		 * Per-genome eggNOG tables ({@code MagDbItem.getOriginalFuncFolder()});
		 * null leaves the {@code assigned_*} function columns NA.
		 */
		public File eggnogDir;
		/** Where the output files go. */
		public File outputDir;
		public int threads = Math.max(1, Math.min(4, Runtime.getRuntime().availableProcessors() - 1));

		/** For the undecidable table logged at the start; the front end used the same values. */
		public double undecidableToleranceDa;
		public double fixedCysteineModDa;

		/** Open search only: MS2 form ranking of the candidates; null leaves them NOT_RANKED. */
		public VariantFormRanking.Inputs ranking;

		public Consumer<String> info = s -> {
		};
		public Consumer<String> warn;
	}

	public static final class Outcome {
		public final boolean ran;
		public final String message;
		public final File candidatesFile;
		public final File sitesFile;
		public final File summaryFile;
		public final File speciesFile;
		public final File columnsFile;

		Outcome(boolean ran, String message, File candidatesFile, File sitesFile, File summaryFile, File speciesFile,
				File columnsFile) {
			this.ran = ran;
			this.message = message;
			this.candidatesFile = candidatesFile;
			this.sitesFile = sitesFile;
			this.summaryFile = summaryFile;
			this.speciesFile = speciesFile;
			this.columnsFile = columnsFile;
		}

		Outcome(boolean ran, String message) {
			this(ran, message, null, null, null, null, null);
		}
	}

	private PangenomeVariantCheck() {
	}

	/**
	 * Whether {@link #run} would do any work with these settings - for a
	 * workflow's skip-to-completion gate, which must not skip when the add-on's
	 * outputs are still to be written. Cheap: looks at file existence only.
	 */
	public static boolean wouldRun(PangenomeVariantSettings settings, File magVersionDir) {
		if (settings == null || settings.getMode() == PangenomeVariantSettings.Mode.OFF) {
			return false;
		}
		PangenomeResource resource = PangenomeResource.locate(settings.getResourcePath(), magVersionDir);
		return resource != null
				&& (!resource.needsPrepare() || settings.getMode() == PangenomeVariantSettings.Mode.ON);
	}

	public static Outcome run(Request r) {
		Consumer<String> info = r.info == null ? s -> {
		} : r.info;
		Consumer<String> warn = r.warn == null ? info : r.warn;
		PangenomeVariantSettings settings = r.settings == null ? new PangenomeVariantSettings() : r.settings;
		PangenomeVariantSettings.Mode mode = settings.getMode();
		if (mode == PangenomeVariantSettings.Mode.OFF) {
			info.accept(LOG_PREFIX + "off (" + PangenomeVariantSettings.KEY_MODE + " = off); skipped");
			return new Outcome(false, "off");
		}
		PangenomeResource resource = PangenomeResource.locate(settings.getResourcePath(), r.magVersionDir);
		if (resource == null) {
			String where = settings.getResourcePath().isEmpty()
					? (r.magVersionDir == null ? "(no MAG database selected)"
							: new File(r.magVersionDir, PangenomeResource.DEFAULT_SUBDIR).getPath())
					: settings.getResourcePath();
			String msg = "UHGP-100 resource not found at " + where + "; skipped. Add it with 'Prepare protein "
					+ "catalogue' in the MAG database configuration, or mag-download --protein-catalogue-only";
			if (mode == PangenomeVariantSettings.Mode.ON) {
				warn.accept(LOG_PREFIX + msg + " (" + PangenomeVariantSettings.KEY_MODE + " = on)");
			} else {
				info.accept(LOG_PREFIX + msg);
			}
			return new Outcome(false, msg);
		}
		try {
			if (resource.needsPrepare()) {
				if (mode != PangenomeVariantSettings.Mode.ON) {
					String msg = "found " + resource.archiveToPrepare + " but it has not been prepared yet (one-time, "
							+ "about 45 min); use 'Prepare protein catalogue' in the MAG database configuration, or "
							+ "set " + PangenomeVariantSettings.KEY_MODE + " = on to prepare it in the next run; "
							+ "skipped";
					info.accept(LOG_PREFIX + msg);
					return new Outcome(false, msg);
				}
				info.accept(LOG_PREFIX + "preparing " + resource.archiveToPrepare + " (one-time)");
				resource = resource.prepare(info);
			}
			return execute(r, settings, resource, info, warn);
		} catch (IOException | RuntimeException e) {
			String msg = "failed, workflow results are unaffected: " + e;
			warn.accept(LOG_PREFIX + msg);
			return new Outcome(false, msg);
		}
	}

	private static Outcome execute(Request r, PangenomeVariantSettings settings, PangenomeResource resource,
			Consumer<String> info, Consumer<String> warn) throws IOException {
		long start = System.currentTimeMillis();
		String engine = engineLabel(r.engine);
		info.accept(LOG_PREFIX + engine + ", " + r.candidates.size() + " substitution candidates, UHGP-100 at "
				+ resource.describe());
		logUndecidableTable(r, info);

		List<VariantCandidate> candidates = new ArrayList<>(r.candidates);
		candidates.sort(CANDIDATE_ORDER);
		List<MatchedControlGenerator.Control> controls = MatchedControlGenerator.generate(candidates,
				r.identifiedBasePeptides, settings.getControlSeed(), settings.getControlsPerCandidate());

		Map<String, Integer> queryIndex = new LinkedHashMap<>();
		for (VariantCandidate c : candidates) {
			queryIndex.putIfAbsent(PeptideAutomaton.normalise(c.variantPeptide), queryIndex.size());
			queryIndex.putIfAbsent(PeptideAutomaton.normalise(c.basePeptide), queryIndex.size());
		}
		for (MatchedControlGenerator.Control c : controls) {
			queryIndex.putIfAbsent(PeptideAutomaton.normalise(c.peptide), queryIndex.size());
		}
		List<String> queries = new ArrayList<>(queryIndex.keySet());
		info.accept(String.format(Locale.ROOT, LOG_PREFIX + "%,d controls (seed %d, %d per candidate); %,d distinct "
				+ "query peptides; tryptic context %s", controls.size(), settings.getControlSeed(),
				settings.getControlsPerCandidate(), queries.size(),
				settings.isRequireTrypticContext() ? "required" : "not required"));

		List<File> searched = new ArrayList<>();
		for (File f : r.searchedFastas) {
			if (f != null && f.isFile()) {
				searched.add(f);
			}
		}
		PangenomeVariantLookup.QueryHit[] searchedHits = null;
		if (searched.isEmpty()) {
			warn.accept(LOG_PREFIX + "no searched database found; in_searched_db is NA");
		} else {
			Map<Integer, java.util.Set<String>> locateIn = new java.util.HashMap<>();
			for (VariantCandidate c : candidates) {
				if (!c.assignedProteins.isEmpty()) {
					locateIn.computeIfAbsent(queryIndex.get(PeptideAutomaton.normalise(c.basePeptide)),
							k -> new java.util.HashSet<>()).addAll(c.assignedProteins);
				}
			}
			searchedHits = PangenomeVariantLookup.scan(queries, searched, null, settings.isRequireTrypticContext(),
					r.threads, r.decoyPrefix, info, "searched database (" + searched.size() + " file(s))", null,
					locateIn);
			locate(candidates, queryIndex, searchedHits, info);
		}

		SpeciesMap species = null;
		if (resource.metadataTsv != null) {
			species = SpeciesMap.load(resource.metadataTsv, resource.repSpeciesGz, info);
		} else {
			warn.accept(LOG_PREFIX + "no " + PangenomeResource.METADATA + " found; n_species is NA and only the "
					+ "ALL scope is reported");
		}
		BitSet selected = species == null ? null : selectedSpecies(r.taxonomyDir, species, info, warn);
		BitSet detected = settings.getSpeciesScope().reportsDetected() ? selected : null;
		File membershipFile = UhgpMembership.locate(settings.getMembershipPath(), resource.directory);
		if (membershipFile == null) {
			String msg = "no " + UhgpMembership.FILE + " cluster membership "
					+ (settings.getMembershipPath().isEmpty() ? "beside " + resource.faaGz
							: "at " + settings.getMembershipPath())
					+ "; genome counts are NA";
			if (settings.getMembershipPath().isEmpty()) {
				info.accept(LOG_PREFIX + msg);
			} else {
				warn.accept(LOG_PREFIX + msg);
			}
		}
		BitSet variantQueries = null;
		if (membershipFile != null) {
			variantQueries = new BitSet();
			for (VariantCandidate c : candidates) {
				variantQueries.set(queryIndex.get(PeptideAutomaton.normalise(c.variantPeptide)));
			}
		}
		PangenomeVariantLookup.QueryHit[] uhgp = PangenomeVariantLookup.scan(queries,
				List.of(resource.faaGz), species, settings.isRequireTrypticContext(), r.threads, "", info,
				"UHGP-100", variantQueries);
		UhgpMembership membership = null;
		if (membershipFile != null) {
			membership = readMembership(membershipFile, uhgp, variantQueries, info, warn);
		}

		if (r.engine == VariantCandidate.Engine.OPEN_SEARCH) {
			VariantFormRanking.run(candidates, r.ranking, r.outputDir, info, warn);
		}

		Map<String, Integer> perPsm = new java.util.HashMap<>();
		for (VariantCandidate c : candidates) {
			perPsm.merge(psmKey(c), 1, Integer::sum);
		}
		for (VariantCandidate c : candidates) {
			c.candidatesForPsm = perPsm.get(psmKey(c));
		}

		Evidence ev = new Evidence(queryIndex, searchedHits, uhgp, detected, settings.getSpeciesScope());
		Map<Integer, int[]> carrierCache = new java.util.HashMap<>();
		AssignedProteinFunction function = AssignedProteinFunction.load(r.eggnogDir, candidates, info, warn);
		for (VariantCandidate c : candidates) {
			if (function != null) {
				function.annotate(c);
			}
			int q = queryIndex.get(PeptideAutomaton.normalise(c.variantPeptide));
			c.inSearchedDb = ev.inSearched(q);
			c.inPangenomeAll = uhgp[q].found();
			c.inPangenomeDetected = ev.inDetected(q);
			c.nSpecies = uhgp[q].species == null ? -1 : uhgp[q].species.cardinality();
			c.exampleProtein = uhgp[q].exampleProtein;
			c.inSameSpecies = sameSpecies(c, species, uhgp[q]);
			if (membership != null) {
				countGenomes(c, species, carriers(q, uhgp[q], membership, carrierCache));
			}
			c.classify();
		}
		logSameSpecies(candidates, info);
		if (membership != null) {
			logGenomeCounts(candidates, info);
		}

		r.outputDir.mkdirs();
		File columnsFile = new File(r.outputDir, COLUMNS_FILE);
		writeAtomically(columnsFile, VariantColumnDictionary::write);
		File candidatesFile = new File(r.outputDir, CANDIDATES_FILE);
		File sitesFile = new File(r.outputDir, SITES_FILE);
		File summaryFile = new File(r.outputDir, SUMMARY_FILE);
		writeAtomically(candidatesFile, w -> writeCandidates(w, candidates));
		List<VariantSiteTable.Row> siteRows = VariantSiteTable.summarise(candidates);
		SpeciesMap siteSpecies = species;
		boolean ranked = r.engine == VariantCandidate.Engine.OPEN_SEARCH;
		writeAtomically(sitesFile, w -> VariantSiteTable.write(w, siteRows, siteSpecies, ranked));
		List<VariantEvidenceSummary.Row> rows = VariantEvidenceSummary.summarise(candidates, controls, ev,
				r.engine == VariantCandidate.Engine.OPEN_SEARCH);
		writeAtomically(summaryFile, w -> VariantEvidenceSummary.write(w, rows));
		File speciesFile = new File(r.outputDir, SPECIES_FILE);
		UhgpMembership members = membership;
		List<VariantSpeciesSummary.Row> speciesRows = species == null ? List.of()
				: VariantSpeciesSummary.summarise(candidates, species, selected, ev,
						members == null ? null : q -> carriers(q, uhgp[q], members, carrierCache));
		SpeciesMap speciesMap = species;
		writeAtomically(speciesFile, w -> VariantSpeciesSummary.write(w, engine, speciesRows, speciesMap));
		if (species == null) {
			warn.accept(LOG_PREFIX + "no genome metadata, " + SPECIES_FILE + " has a header only");
		}

		Map<VariantCandidate.PangenomeClass, Integer> counts = new LinkedHashMap<>();
		for (VariantCandidate.PangenomeClass k : VariantCandidate.PangenomeClass.values()) {
			counts.put(k, 0);
		}
		for (VariantCandidate c : candidates) {
			counts.merge(c.pangenomeClass, 1, Integer::sum);
		}
		info.accept(String.format(Locale.ROOT, LOG_PREFIX + "done in %.0f s; %s; wrote %s, %s, %s, %s and %s",
				(System.currentTimeMillis() - start) / 1000.0, counts, candidatesFile, sitesFile, summaryFile,
				speciesFile, columnsFile));
		return new Outcome(true, "ok", candidatesFile, sitesFile, summaryFile, speciesFile, columnsFile);
	}

	/** One PSM: the candidates of one spectrum's identification share it. */
	private static String psmKey(VariantCandidate c) {
		return String.join("\u0001", String.valueOf(c.file), String.valueOf(c.scan), String.valueOf(c.psmId),
				c.basePeptide);
	}

	/** Copies each candidate's base-peptide locations in its own assigned proteins. */
	private static void locate(List<VariantCandidate> candidates, Map<String, Integer> queryIndex,
			PangenomeVariantLookup.QueryHit[] searchedHits, Consumer<String> info) {
		int pairs = 0;
		int located = 0;
		for (VariantCandidate c : candidates) {
			c.proteinLocations.clear();
			Map<String, PangenomeVariantLookup.Location> at = searchedHits[queryIndex
					.get(PeptideAutomaton.normalise(c.basePeptide))].locations;
			for (String protein : c.assignedProteins) {
				pairs++;
				PangenomeVariantLookup.Location l = at == null ? null : at.get(protein);
				if (l != null) {
					c.proteinLocations.put(protein, l);
					located++;
				}
			}
		}
		info.accept(String.format(Locale.ROOT, LOG_PREFIX + "base peptide located in %,d of %,d (candidate, assigned "
				+ "protein) pairs in the searched databases", located, pairs));
	}

	/** Lookup answers per normalised peptide, shared by candidates, summary and controls. */
	static final class Evidence {
		final Map<String, Integer> queryIndex;
		final PangenomeVariantLookup.QueryHit[] searched;
		final PangenomeVariantLookup.QueryHit[] uhgp;
		final BitSet detectedSpecies;
		final PangenomeVariantSettings.SpeciesScope scope;

		Evidence(Map<String, Integer> queryIndex, PangenomeVariantLookup.QueryHit[] searched,
				PangenomeVariantLookup.QueryHit[] uhgp, BitSet detectedSpecies,
				PangenomeVariantSettings.SpeciesScope scope) {
			this.queryIndex = queryIndex;
			this.searched = searched;
			this.uhgp = uhgp;
			this.detectedSpecies = detectedSpecies;
			this.scope = scope;
		}

		Boolean inSearched(int q) {
			return searched == null ? null : searched[q].found();
		}

		boolean inAll(int q) {
			return uhgp[q].found();
		}

		/** Null when the DETECTED scope is not reported or not available. */
		Boolean inDetected(int q) {
			if (detectedSpecies == null || uhgp[q].species == null) {
				return null;
			}
			return uhgp[q].species.intersects(detectedSpecies);
		}

		boolean reportsAll() {
			return scope.reportsAll();
		}

		boolean reportsDetected() {
			return detectedSpecies != null;
		}

		int index(String peptide) {
			return queryIndex.get(PeptideAutomaton.normalise(peptide));
		}
	}

	/**
	 * Fills {@link VariantCandidate#assignedSpecies} and answers whether any
	 * UHGP-100 hit of the variant spans one of them. Null when there is no
	 * species map, no assigned protein resolves to a species, or the variant's
	 * hits carry no species.
	 */
	static Boolean sameSpecies(VariantCandidate c, SpeciesMap species, PangenomeVariantLookup.QueryHit hit) {
		c.assignedSpecies.clear();
		if (species == null) {
			return null;
		}
		BitSet own = ownSpecies(c, species);
		for (int s = own.nextSetBit(0); s >= 0; s = own.nextSetBit(s + 1)) {
			c.assignedSpecies.add(species.speciesId(s));
		}
		if (own.isEmpty() || hit.species == null) {
			return null;
		}
		return hit.species.intersects(own);
	}

	/** Species indices of the candidate's assigned proteins. */
	static BitSet ownSpecies(VariantCandidate c, SpeciesMap species) {
		BitSet own = new BitSet();
		for (String protein : c.assignedProteins) {
			int s = species.speciesOfAccession(protein);
			if (s >= 0) {
				own.set(s);
			}
		}
		return own;
	}

	/**
	 * Member genomes of every cluster that holds a variant peptide. Null when
	 * the table cannot be read; the genome counts are then NA.
	 */
	private static UhgpMembership readMembership(File file, PangenomeVariantLookup.QueryHit[] uhgp,
			BitSet variantQueries, Consumer<String> info, Consumer<String> warn) {
		int n = 0;
		for (int q = variantQueries.nextSetBit(0); q >= 0; q = variantQueries.nextSetBit(q + 1)) {
			n += uhgp[q].proteinKeyCount;
		}
		long[] keys = new long[n];
		int at = 0;
		for (int q = variantQueries.nextSetBit(0); q >= 0; q = variantQueries.nextSetBit(q + 1)) {
			System.arraycopy(uhgp[q].proteinKeys == null ? new long[0] : uhgp[q].proteinKeys, 0, keys, at,
					uhgp[q].proteinKeyCount);
			at += uhgp[q].proteinKeyCount;
		}
		try {
			UhgpMembership m = UhgpMembership.read(file, keys, info);
			if (m.representativesFound < m.representatives.length) {
				warn.accept(String.format(Locale.ROOT, LOG_PREFIX + "%,d of %,d clusters holding a variant are not "
						+ "in %s; is it from the same catalogue version as %s? Their genomes are not counted",
						m.representatives.length - m.representativesFound, m.representatives.length, file,
						PangenomeResource.FAA));
			}
			return m;
		} catch (IOException e) {
			warn.accept(LOG_PREFIX + "cannot read " + file + ", genome counts are NA: " + e.getMessage());
			return null;
		}
	}

	/** Distinct genomes, sorted, over the clusters a query hit. */
	private static int[] carriers(int q, PangenomeVariantLookup.QueryHit hit, UhgpMembership membership,
			Map<Integer, int[]> cache) {
		return cache.computeIfAbsent(q, k -> {
			int[] all = new int[0];
			int n = 0;
			for (int i = 0; i < hit.proteinKeyCount; i++) {
				int[] g = membership.genomesOf(hit.proteinKeys[i]);
				if (g == null || g.length == 0) {
					continue;
				}
				if (n + g.length > all.length) {
					all = java.util.Arrays.copyOf(all, Math.max(n + g.length, all.length * 2));
				}
				System.arraycopy(g, 0, all, n, g.length);
				n += g.length;
			}
			return java.util.Arrays.stream(all, 0, n).sorted().distinct().toArray();
		});
	}

	/** Sets the genome-count fields from the variant's carrier genomes. */
	static void countGenomes(VariantCandidate c, SpeciesMap species, int[] carriers) {
		c.nGenomes = carriers.length;
		if (species == null) {
			return;
		}
		BitSet own = ownSpecies(c, species);
		if (own.isEmpty()) {
			return;
		}
		int inSpecies = 0;
		for (int s = own.nextSetBit(0); s >= 0; s = own.nextSetBit(s + 1)) {
			inSpecies += species.genomesInSpecies(s);
		}
		int same = 0;
		int isolates = 0;
		long firstGenome = -1;
		long firstIsolate = -1;
		for (int g : carriers) {
			int s = species.speciesOfGenomeNumber(g);
			if (s < 0 || !own.get(s)) {
				continue;
			}
			same++;
			if (firstGenome < 0) {
				firstGenome = g;
			}
			if (species.isIsolate(g)) {
				isolates++;
				if (firstIsolate < 0) {
					firstIsolate = g;
				}
			}
		}
		c.nGenomesInSpecies = inSpecies;
		c.nGenomesSameSpecies = same;
		c.nIsolatesSameSpecies = species.hasGenomeType() ? isolates : -1;
		long example = firstIsolate >= 0 ? firstIsolate : firstGenome;
		c.exampleGenomeSameSpecies = example < 0 ? "" : SpeciesMap.genomeIdOf(example);
	}

	private static void logGenomeCounts(List<VariantCandidate> candidates, Consumer<String> info) {
		int carried = 0;
		int withIsolate = 0;
		for (VariantCandidate c : candidates) {
			if (c.nGenomesSameSpecies > 0) {
				carried++;
				if (c.nIsolatesSameSpecies > 0) {
					withIsolate++;
				}
			}
		}
		info.accept(String.format(Locale.ROOT, LOG_PREFIX + "%,d candidates' variant is carried by at least one "
				+ "genome of their own species, %,d of them by an isolate", carried, withIsolate));
	}

	private static void logSameSpecies(List<VariantCandidate> candidates, Consumer<String> info) {
		int same = 0;
		int other = 0;
		int unknown = 0;
		for (VariantCandidate c : candidates) {
			if (c.inSameSpecies == null) {
				unknown++;
			} else if (c.inSameSpecies) {
				same++;
			} else if (c.inPangenomeAll) {
				other++;
			}
		}
		info.accept(String.format(Locale.ROOT, LOG_PREFIX + "%,d candidates' variant is in UHGP-100 within the "
				+ "species of their own protein, %,d only in other species, %,d with the assigned species unknown",
				same, other, unknown));
	}

	static final Comparator<VariantCandidate> CANDIDATE_ORDER = Comparator
			.comparing((VariantCandidate c) -> c.file == null ? "" : c.file).thenComparingInt(c -> c.scan)
			.thenComparingInt(c -> c.psmId).thenComparing(c -> c.basePeptide).thenComparingInt(c -> c.position)
			.thenComparing(c -> String.valueOf(c.toAa));

	/**
	 * Species representatives of the genomes the run's post-analysis selected
	 * ({@code genome_final_selected.tsv}): the DETECTED scope and the species
	 * summary's {@code in_genome_final_selected} column. Null when the file is missing or names
	 * no known genome; both are then not reported.
	 */
	static BitSet selectedSpecies(File taxonomyDir, SpeciesMap species, Consumer<String> info,
			Consumer<String> warn) throws IOException {
		File f = taxonomyDir == null ? null : new File(taxonomyDir, "genome_final_selected.tsv");
		if (f == null || !f.isFile()) {
			warn.accept(LOG_PREFIX + "no genome_final_selected.tsv in " + taxonomyDir
					+ "; DETECTED scope and selected species not reported");
			return null;
		}
		BitSet out = new BitSet();
		int genomes = 0;
		int unknown = 0;
		try (BufferedReader br = new BufferedReader(new FileReader(f, StandardCharsets.UTF_8))) {
			String line;
			boolean header = true;
			while ((line = br.readLine()) != null) {
				if (line.isBlank()) {
					break;
				}
				if (line.startsWith("#")) {
					continue;
				}
				if (header) {
					header = false;
					if (line.startsWith("Genome\t") || line.equals("Genome")) {
						continue;
					}
				}
				int tab = line.indexOf('\t');
				String genome = (tab < 0 ? line : line.substring(0, tab)).trim();
				genomes++;
				int s = species.speciesOfGenome(genome);
				if (s < 0) {
					unknown++;
				} else {
					out.set(s);
				}
			}
		}
		if (out.isEmpty()) {
			warn.accept(LOG_PREFIX + f + " names no genome of the UHGG metadata; DETECTED scope and selected "
					+ "species not reported");
			return null;
		}
		info.accept(String.format(Locale.ROOT, LOG_PREFIX + "selected species = %d from %d selected genomes%s",
				out.cardinality(), genomes, unknown > 0 ? " (" + unknown + " not in the UHGG metadata)" : ""));
		return out;
	}

	private static void logUndecidableTable(Request r, Consumer<String> info) {
		List<String> table = IsobaricModificationTable.describe(r.undecidableToleranceDa, r.fixedCysteineModDa);
		info.accept(String.format(Locale.ROOT, LOG_PREFIX + "%d substitutions are MS2-undecidable (same-residue "
				+ "modification within %.4f Da, C %+.5f):", table.size(), r.undecidableToleranceDa,
				r.fixedCysteineModDa));
		for (String line : table) {
			info.accept("  " + line);
		}
	}

	static String engineLabel(VariantCandidate.Engine engine) {
		return engine == null ? NA : engine.name().toLowerCase(Locale.ROOT);
	}

	static final String[] CANDIDATE_COLUMNS = { "engine", "file", "scan", "hit_rank", "charge", "precursor_mz",
			"psm_q_value", "base_peptide", "position", "from_aa", "to_aa", "variant_peptide",
			"n_substitutions_in_variant_peptide", "n_candidates_for_psm", "delta_mass", "assigned_proteins",
			"protein_sites", "residue_before", "residue_after", "peak_id", "peak_substitution_shaped",
			"ms2_decidable", "isobaric_modifications", "in_searched_db", "in_pangenome_all",
			"in_pangenome_detected_species", "in_same_species", "n_species", "assigned_species", "example_protein",
			"n_genomes", "n_genomes_same_species", "n_genomes_in_species", "carrier_fraction_in_species",
			"n_isolates_same_species", "example_genome_same_species", "ms2_verdict", "ms2_score_substitution",
			"ms2_score_other_substitution", "ms2_score_modification", "ms2_score_elsewhere",
			"ms2_localises_elsewhere", "ms2_best_form", "pangenome_class", "assigned_kegg_ko", "assigned_cog",
			"assigned_cog_category", "assigned_preferred_name", "assigned_description", "front_end_call",
			"localized_sequence", "reassigner_in_db" };

	static void writeCandidates(Writer w, List<VariantCandidate> candidates) throws IOException {
		w.write(String.join("\t", CANDIDATE_COLUMNS));
		w.write('\n');
		for (VariantCandidate c : candidates) {
			boolean openSearch = c.engine == VariantCandidate.Engine.OPEN_SEARCH;
			List<String> sites = new ArrayList<>();
			PangenomeVariantLookup.Location flank = null;
			for (Map.Entry<String, PangenomeVariantLookup.Location> e : c.proteinLocations.entrySet()) {
				sites.add(e.getKey() + ":" + c.proteinPosition(e.getKey()));
				if (flank == null) {
					flank = e.getValue();
				}
			}
			int nSubstitutions = c.substitutionsInVariantPeptide();
			String[] f = { engineLabel(c.engine), clean(c.file), String.valueOf(c.scan),
					openSearch ? String.valueOf(c.psmId) : NA, c.charge > 0 ? String.valueOf(c.charge) : NA,
					Double.isNaN(c.precursorMz) ? NA : String.format(Locale.ROOT, "%.5f", c.precursorMz),
					Double.isNaN(c.qValue) ? NA : String.format(Locale.ROOT, "%.3e", c.qValue), c.basePeptide,
					String.valueOf(c.position), String.valueOf(c.fromAa), String.valueOf(c.toAa), c.variantPeptide,
					nSubstitutions < 0 ? NA : String.valueOf(nSubstitutions), String.valueOf(c.candidatesForPsm),
					String.format(Locale.ROOT, "%.5f", c.deltaMass),
					c.assignedProteins.isEmpty() ? NA : clean(String.join(";", c.assignedProteins)),
					sites.isEmpty() ? NA : clean(String.join(";", sites)),
					flank == null ? NA : String.valueOf(flank.before), flank == null ? NA : String.valueOf(flank.after),
					c.peakId < 0 ? NA : String.valueOf(c.peakId), bool(c.peakSubstitutionShaped),
					String.valueOf(c.ms2Decidable), c.isobaricModifications.isEmpty() ? ""
							: clean(String.join("; ", c.isobaricModifications)),
					bool(c.inSearchedDb), String.valueOf(c.inPangenomeAll), bool(c.inPangenomeDetected),
					bool(c.inSameSpecies), c.nSpecies < 0 ? NA : String.valueOf(c.nSpecies),
					c.assignedSpecies.isEmpty() ? NA : String.join(";", c.assignedSpecies), clean(c.exampleProtein),
					count(c.nGenomes), count(c.nGenomesSameSpecies), count(c.nGenomesInSpecies), carrierFraction(c),
					count(c.nIsolatesSameSpecies), c.exampleGenomeSameSpecies.isEmpty() ? NA : c.exampleGenomeSameSpecies,
					c.ms2Verdict.name(),
					score(c.ms2ScoreSubstitution), score(c.ms2ScoreOtherSubstitution), score(c.ms2ScoreModification),
					score(c.ms2ScoreElsewhere), bool(c.ms2LocalisesElsewhere),
					c.ms2BestForm.isEmpty() ? NA : clean(c.ms2BestForm),
					c.pangenomeClass == null ? NA : c.pangenomeClass.name(), function(c.assignedKeggKo),
					function(c.assignedCog), function(c.assignedCogCategory), function(c.assignedPreferredName),
					function(c.assignedDescription), c.frontEndEvidence.isEmpty() ? NA : clean(c.frontEndEvidence),
					openSearch && !c.localizedSequence.isEmpty() ? clean(c.localizedSequence) : NA,
					bool(c.reassignerInDb) };
			w.write(String.join("\t", f));
			w.write('\n');
		}
	}

	/** {@code n_genomes_same_species / n_genomes_in_species}; NA when either is unknown or the species is empty. */
	static String carrierFraction(VariantCandidate c) {
		if (c.nGenomesSameSpecies < 0 || c.nGenomesInSpecies <= 0) {
			return NA;
		}
		return String.format(Locale.ROOT, "%.4f", (double) c.nGenomesSameSpecies / c.nGenomesInSpecies);
	}

	static String bool(Boolean b) {
		return b == null ? NA : b.toString();
	}

	private static String score(double s) {
		return Double.isNaN(s) ? NA : String.format(Locale.ROOT, "%.4f", s);
	}

	private static String count(int n) {
		return n < 0 ? NA : String.valueOf(n);
	}

	private static String function(String s) {
		return s == null ? NA : clean(s);
	}

	static String clean(String s) {
		return s == null ? "" : s.replace('\t', ' ').replace('\n', ' ').replace('\r', ' ');
	}

	interface WriterBody {
		void write(Writer w) throws IOException;
	}

	/** Writes beside the target and renames, so a crash never leaves a truncated file that looks complete. */
	static void writeAtomically(File target, WriterBody body) throws IOException {
		File part = new File(target.getParentFile(), target.getName() + ".part");
		try (Writer w = new BufferedWriter(
				new OutputStreamWriter(Files.newOutputStream(part.toPath()), StandardCharsets.UTF_8))) {
			body.write(w);
		}
		Files.move(part.toPath(), target.toPath(), StandardCopyOption.REPLACE_EXISTING);
	}
}
