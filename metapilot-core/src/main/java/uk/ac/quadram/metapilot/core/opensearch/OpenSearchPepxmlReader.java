package uk.ac.quadram.metapilot.core.opensearch;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.Set;
import java.util.function.Consumer;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.dom4j.Attribute;
import org.dom4j.DocumentException;
import org.dom4j.ElementHandler;
import org.dom4j.ElementPath;
import org.dom4j.Element;
import org.dom4j.io.SAXReader;

/**
 * Reads MSFragger's native pepXML for the standalone open-search pipeline.
 * <p>
 * Differs from {@code dbSearch.msFragger.MSFraggerReader} in three ways that
 * are each required here rather than stylistic:
 * <ol>
 * <li><b>Every candidate, not just the top hit.</b> The existing reader takes
 * only the first {@code search_hit} per {@code spectrum_query}. With
 * {@code output_report_topN > 1} that silently discards ranks 2..N, which
 * would defeat this pipeline's requirement that the classifier decide, not the
 * engine.</li>
	 * <li><b>Localization is read from {@code ptm_result}.</b> It is a sibling of
	 * {@code search_score}, not one of the search-score entries. Both MSFragger's
	 * legacy compact attributes and its expanded 4.1+ attributes are mapped.</li>
	 * <li><b>The emitted vocabulary is retained and reported.</b> Every search
	 * score and every {@code ptm_result} attribute is captured before known fields
	 * are mapped, so a future MSFragger rename is visible in the run log rather
	 * than silently turning classifier features into NaN.</li>
 * </ol>
 */
public final class OpenSearchPepxmlReader {

	private static final Logger LOGGER = LogManager.getLogger(OpenSearchPepxmlReader.class);

	/**
	 * Fallback {@code search_score} aliases retained for compatibility with
	 * converted/non-native pepXML. Native MSFragger output uses
	 * {@code ptm_result}; see {@link #readPtmResult(Element, OpenSearchPsm)}.
	 */
	private static final List<String> LOCALIZED_SEQUENCE_ALIASES = Arrays.asList("best_locs", "msfragger localization",
			"best_positions", "best positions");
	private static final List<String> SCORE_WITHOUT_DELTA_ALIASES = Arrays.asList("score_without_delta_mass",
			"score all unshifted", "score_all_unshifted");
	private static final List<String> SCORE_WITH_DELTA_ALIASES = Arrays.asList("best_score_with_delta_mass",
			"score shifted best position", "score_shifted_best_position", "score best position",
			"score_best_position");
	private static final List<String> SECOND_SCORE_WITH_DELTA_ALIASES = Arrays.asList(
			"second_best_score_with_delta_mass", "score second best position", "score_second_best_position");
	private static final List<String> IONS_WITHOUT_DELTA_ALIASES = Arrays.asList("ions all unshifted",
			"ions_all_unshifted");
	private static final List<String> IONS_WITH_DELTA_ALIASES = Arrays.asList("ions shifted best position",
			"ions_shifted_best_position", "ions best position", "ions_best_position");

	/** Score names already handled; anything else is reported once as unknown. */
	private static final Set<String> KNOWN_SCORES = buildKnownScores();

	/** How often a long read reports progress, in spectra. */
	private static final int PROGRESS_EVERY_SPECTRA = 20000;

	private final String decoyTag;

	/**
	 * Highest {@code hit_rank} to retain; 0 or less means keep every rank.
	 * <p>
	 * Stage 1 sets this to 1: it maps genomes from unmodified peptides and never
	 * arbitrates between competing candidates, so ranks 2..N are pure overhead
	 * there. Stage 2 leaves it unset, because adjudicating those candidates is
	 * precisely the classifier's job.
	 */
	private int maxHitRank = 0;

	public void setMaxHitRank(int maxHitRank) {
		this.maxHitRank = maxHitRank;
	}

	/**
	 * Reads a closed search: each PSM is marked closed-stage and carries the variable
	 * modifications of its {@code modification_info}. An open search's delta mass is not a
	 * variable modification, so this is off for it.
	 */
	private boolean closedSearch = false;

	public void setClosedSearch(boolean closedSearch) {
		this.closedSearch = closedSearch;
	}
	private final Set<String> unrecognisedScoreNames = new LinkedHashSet<>();

	/**
	 * Every distinct {@code search_score} name encountered, recognised or not.
	 * <p>
	 * Reported unconditionally rather than only on failure. Reconstructing the
	 * vocabulary afterwards would otherwise mean grepping a multi-gigabyte
	 * pepXML.
	 */
	private final Set<String> allScoreNames = new LinkedHashSet<>();
	private final Set<String> allPtmResultAttributeNames = new LinkedHashSet<>();

	public OpenSearchPepxmlReader(String decoyTag) {
		this.decoyTag = decoyTag == null || decoyTag.isBlank() ? "rev_" : decoyTag;
	}

	private static Set<String> buildKnownScores() {
		Set<String> known = new LinkedHashSet<>(Arrays.asList("hyperscore", "nextscore", "expect", "expectscore"));
		known.addAll(LOCALIZED_SEQUENCE_ALIASES);
		known.addAll(SCORE_WITHOUT_DELTA_ALIASES);
		known.addAll(SCORE_WITH_DELTA_ALIASES);
		known.addAll(SECOND_SCORE_WITH_DELTA_ALIASES);
		known.addAll(IONS_WITHOUT_DELTA_ALIASES);
		known.addAll(IONS_WITH_DELTA_ALIASES);
		return known;
	}

	/** Reads every PSM candidate from every supplied pepXML. */
	public List<OpenSearchPsm> read(List<File> pepxmlFiles) throws DocumentException {
		return read(pepxmlFiles, null);
	}

	/**
	 * Reads every PSM candidate from every supplied pepXML, reporting progress.
	 * This is the largest read of a run — all samples pooled for validation — so
	 * it reports per file as well as within each file.
	 */
	public List<OpenSearchPsm> read(List<File> pepxmlFiles, Consumer<String> progress) throws DocumentException {
		List<OpenSearchPsm> all = new ArrayList<>();
		for (int i = 0; i < pepxmlFiles.size(); i++) {
			File file = pepxmlFiles.get(i);
			if (progress != null) {
				progress.accept("Reading pepXML " + (i + 1) + "/" + pepxmlFiles.size() + ": " + file.getName()
						+ " (" + (file.length() / (1024 * 1024)) + " MB)");
			}
			all.addAll(read(file, progress));
		}
		reportScoreNames(progress);
		return all;
	}

	/**
	 * Streams the pepXML rather than building a whole DOM.
	 * <p>
	 * An open search with {@code output_report_topN > 1} produces a pepXML
	 * several times larger than a closed search's, and loading one into memory
	 * as a document tree costs roughly an order of magnitude more heap than the
	 * file itself — enough to spend minutes in GC or fail outright on a large
	 * sample, with no outward sign of progress. Each {@code spectrum_query} is
	 * therefore parsed and detached as it is encountered, which keeps memory
	 * bounded regardless of file size.
	 *
	 * @param progress optional sink for periodic "N spectra read" updates; a long
	 *                 read is otherwise silent and indistinguishable from a hang
	 */
	public List<OpenSearchPsm> read(File pepxmlFile, Consumer<String> progress) throws DocumentException {
		List<OpenSearchPsm> psms = new ArrayList<>();
		String[] baseName = { baseNameOf(null, pepxmlFile) };
		int[] spectraSeen = { 0 };

		SAXReader reader = new SAXReader();

		// base_name is an attribute of the start tag, so it is available before
		// any spectrum_query beneath it is handled.
		reader.addHandler("/msms_pipeline_analysis/msms_run_summary", new ElementHandler() {
			@Override
			public void onStart(ElementPath path) {
				String raw = path.getCurrent().attributeValue("base_name");
				baseName[0] = baseNameOf(raw, pepxmlFile);
			}

			@Override
			public void onEnd(ElementPath path) {
				path.getCurrent().detach();
			}
		});

		reader.addHandler("/msms_pipeline_analysis/msms_run_summary/spectrum_query", new ElementHandler() {
			@Override
			public void onStart(ElementPath path) {
				// Nothing: the element is not yet populated.
			}

			@Override
			public void onEnd(ElementPath path) {
				Element query = path.getCurrent();
				readSpectrumQuery(query, baseName[0], psms);
				spectraSeen[0]++;
				if (progress != null && spectraSeen[0] % PROGRESS_EVERY_SPECTRA == 0) {
					progress.accept(pepxmlFile.getName() + ": " + spectraSeen[0] + " spectra read, " + psms.size()
							+ " PSM candidates so far ...");
				}
				// Detaching is what bounds memory - without it the handler still
				// accumulates the full tree.
				query.detach();
			}
		});

		reader.read(pepxmlFile);

		LOGGER.info("Read {} PSM candidates from {} spectra in {}", psms.size(), spectraSeen[0],
				pepxmlFile.getName());
		return psms;
	}

	public List<OpenSearchPsm> read(File pepxmlFile) throws DocumentException {
		return read(pepxmlFile, null);
	}

	private void readSpectrumQuery(Element query, String baseName, List<OpenSearchPsm> out) {
		int scan = parseInt(query.attributeValue("start_scan"), -1);
		int charge = parseInt(query.attributeValue("assumed_charge"), 0);
		double precursorMass = parseDouble(query.attributeValue("precursor_neutral_mass"), Double.NaN);
		double rt = parseDouble(query.attributeValue("retention_time_sec"), Double.NaN);

		Element result = query.element("search_result");
		if (result == null) {
			return;
		}

		// Every search_hit, not just the first: ranks 2..N are exactly the
		// candidates the classifier is meant to adjudicate.
		Iterator<Element> hitIt = result.elementIterator("search_hit");
		while (hitIt.hasNext()) {
			Element hit = hitIt.next();
			// Skipped before the PSM is built, so the lower ranks cost nothing in
			// allocation when a caller (stage 1) only wants the top hit.
			if (maxHitRank > 0 && parseInt(hit.attributeValue("hit_rank"), 1) > maxHitRank) {
				continue;
			}
			OpenSearchPsm psm = new OpenSearchPsm();
			psm.setFileName(baseName);
			psm.setScan(scan);
			psm.setCharge(charge);
			psm.setPrecursorNeutralMass(precursorMass);
			psm.setRetentionTime(rt);

			psm.setHitRank(parseInt(hit.attributeValue("hit_rank"), 1));
			psm.setPeptide(hit.attributeValue("peptide"));
			psm.setPrevAA(hit.attributeValue("peptide_prev_aa"));
			psm.setNextAA(hit.attributeValue("peptide_next_aa"));
			double calcMass = parseDouble(hit.attributeValue("calc_neutral_pep_mass"), Double.NaN);
			psm.setCalcNeutralPepMass(calcMass);
			// Derived, not read. massdiff is by definition precursor_neutral_mass
			// minus calc_neutral_pep_mass, and deriving it is what lets Crystal-C's
			// refinement reach every rank: Crystal-C rewrites the precursor mass on
			// the spectrum_query - the correction is a property of the spectrum,
			// not of a hit - and FragPipe only ever feeds it a top-1 pepXML, so
			// nothing upstream guarantees it rewrote the massdiff attribute on
			// ranks 2..N as well. The identity was checked on both kinds of file:
			// an uncorrected MetaPilot pepXML (659,521 hits, max deviation
			// 4.7e-4 Da, float rounding) and a Crystal-C output (106,340 hits, max
			// deviation 1e-5 Da). So this is a no-op when the step is off and the
			// fix when it is on. The attribute is the fallback for a file that
			// carries no precursor mass at all.
			double derived = precursorMass - calcMass;
			psm.setMassDiff(Double.isNaN(derived) ? parseDouble(hit.attributeValue("massdiff"), Double.NaN)
					: derived);
			psm.setNumTolTerm(parseInt(hit.attributeValue("num_tol_term"), -1));
			psm.setNumMissedCleavages(parseInt(hit.attributeValue("num_missed_cleavages"), -1));
			psm.setNumMatchedIons(parseInt(hit.attributeValue("num_matched_ions"), -1));
			psm.setTotNumIons(parseInt(hit.attributeValue("tot_num_ions"), -1));

			String protein = hit.attributeValue("protein");
			psm.setProtein(protein);
			psm.setAllProteins(readAllProteins(hit, protein));
			// Target/decoy is decided on the delegate only, matching how
			// Philosopher sets "Is Decoy" in psm.tsv. A target peptide can carry
			// decoy entries among its alternatives when sequences coincide;
			// treating those as decoy hits would corrupt the 1:1 FDR.
			psm.setTarget(protein == null || !protein.toLowerCase(Locale.ROOT)
					.startsWith(decoyTag.toLowerCase(Locale.ROOT)));

			readSearchScores(hit, psm);
			readPtmResult(hit, psm);
			if (closedSearch) {
				psm.setClosedStage(true);
				psm.setClosedModifications(readVariableModifications(hit, psm.getPeptide()));
			}
			out.add(psm);
		}
	}

	/** Mass of the free peptide N-terminus' hydrogen, which {@code mod_nterm_mass} includes. */
	private static final double NTERM_HYDROGEN = 1.007825;
	/** Below this a position's mass is taken as unmodified (rounding in the pepXML). */
	private static final double MODIFICATION_EPSILON = 0.01;

	/**
	 * Variable modifications from {@code modification_info}, indexed as
	 * {@link OpenSearchPsm#getClosedModifications()}; null when there are none. A residue's
	 * {@code variable} attribute is used when MSFragger writes it; otherwise the variable part is
	 * its mass less the residue with the fixed cysteine modification.
	 */
	static double[] readVariableModifications(Element hit, String peptide) {
		Element info = hit.element("modification_info");
		if (info == null || peptide == null) {
			return null;
		}
		double[] mods = new double[peptide.length() + 1];
		boolean any = false;
		double nterm = parseDouble(info.attributeValue("mod_nterm_mass"), Double.NaN);
		if (Double.isFinite(nterm) && nterm - NTERM_HYDROGEN > MODIFICATION_EPSILON) {
			mods[0] = nterm - NTERM_HYDROGEN;
			any = true;
		}
		Iterator<Element> it = info.elementIterator("mod_aminoacid_mass");
		while (it.hasNext()) {
			Element mod = it.next();
			int position = parseInt(mod.attributeValue("position"), -1);
			if (position < 1 || position > peptide.length()) {
				continue;
			}
			double variable = parseDouble(mod.attributeValue("variable"), Double.NaN);
			if (Double.isNaN(variable)) {
				double mass = parseDouble(mod.attributeValue("mass"), Double.NaN);
				double residue = SubstitutionMassTable.residueMass(peptide.charAt(position - 1));
				variable = Double.isFinite(mass) && residue > 0 ? mass - residue : 0.0;
			}
			if (Math.abs(variable) > MODIFICATION_EPSILON) {
				mods[position] += variable;
				any = true;
			}
		}
		return any ? mods : null;
	}

	/**
	 * Reads MSFragger's localization payload. It is a sibling of
	 * {@code search_score}, not a search score itself; looking only at
	 * {@code search_score} was why every localization feature was NaN in the
	 * first real run. Attribute names changed between the compact pre-4.1 form
	 * and the expanded 4.1+ form, so both are accepted.
	 */
	private void readPtmResult(Element hit, OpenSearchPsm psm) {
		Element ptm = hit.element("ptm_result");
		if (ptm == null) {
			return;
		}

		Map<String, String> attributes = new LinkedHashMap<>();
		Iterator<Attribute> attributeIt = ptm.attributeIterator();
		while (attributeIt.hasNext()) {
			Attribute attribute = attributeIt.next();
			attributes.put(attribute.getName(), attribute.getValue());
			allPtmResultAttributeNames.add(attribute.getName());
		}
		psm.setPtmResultAttributes(attributes);

		String localized = first(attributes, "localization_peptide");
		if (localized != null) {
			psm.setLocalizedSequence(localized);
		}
		double value = parseDouble(first(attributes, "score_all_unshifted", "score_without_ptm"), Double.NaN);
		if (!Double.isNaN(value)) {
			psm.setScoreWithoutDeltaMass(value);
		}
		value = parseDouble(first(attributes, "score_best_position", "best_score_with_ptm"), Double.NaN);
		if (!Double.isNaN(value)) {
			psm.setScoreWithDeltaMass(value);
		}
		value = parseDouble(first(attributes, "score_second_best_position", "second_best_score_with_ptm"),
				Double.NaN);
		if (!Double.isNaN(value)) {
			psm.setSecondBestScoreWithDeltaMass(value);
		}
		value = parseDouble(first(attributes, "ions_all_unshifted"), Double.NaN);
		if (!Double.isNaN(value)) {
			psm.setIonsWithoutDeltaMass(value);
		}
		value = parseDouble(first(attributes, "ions_shifted_best_position", "ions_best_position"), Double.NaN);
		if (!Double.isNaN(value)) {
			psm.setIonsWithDeltaMass(value);
		}
	}

	private static String first(Map<String, String> values, String... names) {
		for (String name : names) {
			String value = values.get(name);
			if (value != null && !value.isBlank()) {
				return value;
			}
		}
		return null;
	}

	/**
	 * Delegate protein plus every {@code <alternative_protein>} child, in
	 * document order and de-duplicated.
	 * <p>
	 * MAG genome selection counts genomes per peptide, so dropping the
	 * alternatives would silently credit a marker peptide shared across many
	 * genomes to only the delegate one and distort the knee cut. Philosopher's
	 * psm.tsv exposes the same information as "Protein" + "Mapped Proteins";
	 * this reader bypasses Philosopher, so it reconstructs the list itself.
	 */
	private static String[] readAllProteins(Element hit, String delegate) {
		Set<String> proteins = new LinkedHashSet<>();
		if (delegate != null && !delegate.trim().isEmpty()) {
			proteins.add(delegate.trim());
		}
		Iterator<Element> altIt = hit.elementIterator("alternative_protein");
		while (altIt.hasNext()) {
			String alt = altIt.next().attributeValue("protein");
			if (alt != null && !alt.trim().isEmpty()) {
				proteins.add(alt.trim());
			}
		}
		return proteins.toArray(new String[0]);
	}

	private void readSearchScores(Element hit, OpenSearchPsm psm) {
		Map<String, String> scores = new LinkedHashMap<>();
		Iterator<Element> scoreIt = hit.elementIterator("search_score");
		while (scoreIt.hasNext()) {
			Element score = scoreIt.next();
			String name = score.attributeValue("name");
			String value = score.attributeValue("value");
			if (name == null) {
				continue;
			}
			// Stored as a string and converted only where a numeric field wants
			// it: at least one localization score is a peptide sequence, and
			// parsing eagerly is what breaks the pre-existing readers.
			scores.put(name, value);
			allScoreNames.add(name);
		}
		psm.setSearchScores(scores);

		for (Map.Entry<String, String> entry : scores.entrySet()) {
			String key = entry.getKey().toLowerCase(Locale.ROOT).trim();
			String value = entry.getValue();
			if (key.equals("hyperscore")) {
				psm.setHyperscore(parseDouble(value, Double.NaN));
			} else if (key.equals("nextscore")) {
				psm.setNextscore(parseDouble(value, Double.NaN));
			} else if (key.equals("expect") || key.equals("expectscore")) {
				psm.setExpect(parseDouble(value, Double.NaN));
			} else if (LOCALIZED_SEQUENCE_ALIASES.contains(key)) {
				psm.setLocalizedSequence(value);
			} else if (SCORE_WITHOUT_DELTA_ALIASES.contains(key)) {
				psm.setScoreWithoutDeltaMass(parseDouble(value, Double.NaN));
			} else if (SCORE_WITH_DELTA_ALIASES.contains(key)) {
				psm.setScoreWithDeltaMass(parseDouble(value, Double.NaN));
			} else if (SECOND_SCORE_WITH_DELTA_ALIASES.contains(key)) {
				psm.setSecondBestScoreWithDeltaMass(parseDouble(value, Double.NaN));
			} else if (IONS_WITHOUT_DELTA_ALIASES.contains(key)) {
				psm.setIonsWithoutDeltaMass(parseDouble(value, Double.NaN));
			} else if (IONS_WITH_DELTA_ALIASES.contains(key)) {
				psm.setIonsWithDeltaMass(parseDouble(value, Double.NaN));
			} else if (!KNOWN_SCORES.contains(key)) {
				unrecognisedScoreNames.add(entry.getKey());
			}
		}
	}

	/**
	 * Logs any {@code search_score} names this reader did not map. Because the
	 * pepXML vocabulary is undocumented, this is how the first real run tells us
	 * what MSFragger actually emits — the values are still on the PSM via
	 * {@code getSearchScores()}, so nothing is lost meanwhile.
	 */
	public void reportUnrecognisedScores() {
		reportScoreNames(null);
	}

	/**
	 * Reports what the pepXML actually contained, through the progress callback
	 * as well as the class logger.
	 * <p>
	 * The callback is the point. Task logs are written by a non-additive
	 * {@code TaskLogger_<id>}, so a bare {@code LOGGER.info} from this class does
	 * not reach the run's {@code logging.log} at all — a diagnostic nobody can
	 * read is worse than none, because its absence looks like evidence.
	 */
	public void reportScoreNames(Consumer<String> progress) {
		String found = allScoreNames.isEmpty() ? "(none)" : String.join(", ", allScoreNames);
		emit(progress, "MSFragger pepXML search_score names found: " + found);
		String ptmFound = allPtmResultAttributeNames.isEmpty() ? "(none)"
				: String.join(", ", allPtmResultAttributeNames);
		emit(progress, "MSFragger pepXML ptm_result attributes found: " + ptmFound);

		if (!unrecognisedScoreNames.isEmpty()) {
			emit(progress, "Of those, names this reader does not map yet: "
					+ String.join(", ", unrecognisedScoreNames)
					+ ". They are retained per-PSM via getSearchScores(); add aliases in OpenSearchPepxmlReader "
					+ "if any correspond to delta-mass localization fields.");
		}
	}

	private static void emit(Consumer<String> progress, String message) {
		LOGGER.info(message);
		if (progress != null) {
			progress.accept(message);
		}
	}

	public Set<String> getUnrecognisedScoreNames() {
		return Collections.unmodifiableSet(unrecognisedScoreNames);
	}

	/**
	 * Verifies that the localization evidence the classifier depends on is
	 * actually present. Called after reading so a misconfigured search
	 * ({@code localize_delta_mass} off, or an MSFragger version using different
	 * attributes) is reported rather than silently training on NaN features.
	 *
	 * @return true when at least one PSM carries a usable delta-mass score gain
	 */
	public static boolean hasLocalizationEvidence(List<OpenSearchPsm> psms) {
		for (OpenSearchPsm psm : psms) {
			if (!Double.isNaN(psm.getDeltaMassScoreGain())
					|| (psm.getLocalizedSequence() != null && !psm.getLocalizedSequence().isBlank())) {
				return true;
			}
		}
		return false;
	}

	private static String baseNameOf(String baseName, File fallback) {
		if (baseName == null || baseName.isEmpty()) {
			String name = fallback.getName();
			int dot = name.indexOf('.');
			return dot > 0 ? name.substring(0, dot) : name;
		}
		int slash = Math.max(baseName.lastIndexOf('\\'), baseName.lastIndexOf('/'));
		return slash >= 0 && slash < baseName.length() - 1 ? baseName.substring(slash + 1) : baseName;
	}

	private static int parseInt(String value, int fallback) {
		if (value == null || value.isBlank()) {
			return fallback;
		}
		try {
			return Integer.parseInt(value.trim());
		} catch (NumberFormatException e) {
			return fallback;
		}
	}

	private static double parseDouble(String value, double fallback) {
		if (value == null || value.isBlank()) {
			return fallback;
		}
		try {
			return Double.parseDouble(value.trim());
		} catch (NumberFormatException e) {
			return fallback;
		}
	}
}
