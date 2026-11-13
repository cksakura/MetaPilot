package uk.ac.quadram.metapilot.core.opensearch;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.StandardCopyOption;
import java.util.Arrays;
import java.util.List;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import uk.ac.quadram.metapilot.core.mag.MagResultLayoutUtil;
import uk.ac.quadram.metapilot.core.crystalc.CrystalCHeadlessRunner;
import uk.ac.quadram.metapilot.core.msbooster.MSBoosterHeadlessRunner;

/**
 * Single source of truth for the open-search workflow's result-folder layout
 * and its resume sentinels.
 * <p>
 * The runner and {@code ResultFolderCleaner} both need to know which files gate
 * skipping already-completed work. Written independently they would drift, and
 * the failure is silent — so both must read them from here.
 *
 * <h2>The resume contract</h2>
 * Each expensive step declares one sentinel. If it exists and is non-empty the
 * step is skipped.
 * <ol>
 * <li><b>The sentinel must be the last file a step writes.</b> A crash midway
 * then leaves it absent and the step correctly reruns, rather than being
 * falsely skipped with partial outputs on disk.</li>
 * <li><b>Every file a downstream consumer reads must be protected</b> — either
 * by being the sentinel itself or by an explicit cleaner keep-rule. The
 * sentinel's presence suppresses regeneration, so deleting a non-sentinel
 * output the step also produced yields a silent partial failure.</li>
 * </ol>
 * Rule 2 is not hypothetical: FragPipe's {@code combined/<expName>/psm.tsv} is
 * not resume-gated, but post-analysis reads PSM-level fields from it on every
 * run. Deleting those files while keeping {@code combined_modified_peptide.tsv}
 * produces a {@code peptides=N, psms=0} result, because the combined search's
 * own sentinel is still present and the step skips.
 */
public final class OpenSearchLayout {

	/** Parameter filename; also the cleaner's workflow-detection key. */
	public static final String PARAMETER_FILE = "mag_opensearch_analysis_parameters.json";

	/**
	 * Engine-named search folder directly under the user's result folder, the
	 * same convention every other workflow uses ({@code fragpipe_MAG_search},
	 * {@code pfind_MAG_search}, {@code sage_MAG_search}):
	 *
	 * <pre>
	 * opensearch_MAG_search/
	 *   search_results/&lt;sample&gt;/marker|sample_db|sample_db_closed/
	 *   combined/
	 *   sequence_variants/
	 *   taxonomy_analysis/
	 *   functional_annotation/
	 *   excluded_evidence/
	 *   report/
	 * </pre>
	 *
	 * Closed and open MSFragger searches of the refined database replace a
	 * second engine, so there is no nested {@code msfragger_only/} or
	 * {@code msfragger_sage/}. {@code MetaBaseParameter.setDbSearchResultFile}
	 * is this folder, where the shared report generator puts {@code report/}.
	 */
	public static final String ENGINE_SEARCH_DIR = "opensearch_MAG_search";

	/**
	 * Pre-flatten layout: pooled validation lived in
	 * {@code opensearch_MAG_search/msfragger_only/} so MSFragger-only and
	 * MSFragger+Sage results could coexist. Kept so leftover folders are not
	 * treated as sample directories, and so {@link #migrateLegacyModeRoot}
	 * can lift that tree to the engine root.
	 */
	public static final String LEGACY_MSFRAGGER_ONLY_DIR = "msfragger_only";
	public static final String LEGACY_MSFRAGGER_SAGE_DIR = "msfragger_sage";

	/** Lower-case form used by the cleaner's directory matching. */
	public static final String ENGINE_SEARCH_DIR_KEY = "opensearch_mag_search";

	public static final String SEARCH_RESULTS_DIR = MagResultLayoutUtil.SEARCH_RESULTS_DIR_NAME;
	/** Round one — genome discovery. Same name the other MAG engines use. */
	public static final String MARKER_SEARCH_DIR = "marker";
	/** Round two — search against the refined per-sample database. */
	public static final String SAMPLE_DB_SEARCH_DIR = "sample_db";
	/** Closed search of the refined database, before the open one. */
	public static final String SAMPLE_DB_CLOSED_SEARCH_DIR = "sample_db_closed";
	/** Pooled, cross-sample validation output. */
	public static final String COMBINED_DIR = "combined";
	/**
	 * Sequence-variant calls and the UHGP-100 check of those substitutions.
	 * Beside {@code combined/}, not inside it: the identification tables stay
	 * where the other MAG workflows put them.
	 */
	public static final String SEQUENCE_VARIANTS_DIR = "sequence_variants";

	/**
	 * MSFragger parameter file, written into the same folder as the pepXML it
	 * produced.
	 * <p>
	 * Every search folder must be self-describing: opening it should show both
	 * the result and the exact parameters that produced it, with no cross
	 * referencing against a folder elsewhere in the tree. This matters most when
	 * a run is resumed — a skipped step's parameters stay pinned next to its
	 * output, so a later parameter change is visible as a mismatch rather than
	 * silently applying to only the steps that reran.
	 */
	public static final String SEARCH_PARAMS_FILE = "fragger_open.params";

	private static final Logger LOGGER = LogManager.getLogger(OpenSearchLayout.class);

	private OpenSearchLayout() {
	}

	/**
	 * The engine search root inside the user's chosen result folder. All other
	 * paths in this class are relative to this, not to the user's folder
	 * directly — matching how the other workflows nest their output.
	 */
	public static File engineSearchRoot(File userResultFolder) {
		return new File(userResultFolder, ENGINE_SEARCH_DIR);
	}

	/**
	 * Delegates to {@link MagResultLayoutUtil} so this workflow inherits the
	 * same legacy-flat-layout migration the other MAG engines get.
	 */
	public static File searchResultsRoot(File resultRoot) {
		migrateLegacyModeRoot(resultRoot);
		return MagResultLayoutUtil.resolveSearchResultsRoot(resultRoot, workflowDirNames());
	}

	/**
	 * This workflow's own folders beside {@code search_results/}. The legacy
	 * migration moves every other unrecognised folder into it, and did move
	 * {@code combined_msfragger_sage/} that way in earlier result folders.
	 */
	static java.util.Set<String> workflowDirNames() {
		java.util.Set<String> names = new java.util.HashSet<>(modeFolderNames());
		names.add(SEQUENCE_VARIANTS_DIR);
		names.add(MagResultLayoutUtil.SPECTRA_SUMMARY_DIR_NAME);
		names.add(OpenSearchSpectra.CONVERTED_DIR);
		return names;
	}

	public static File sampleDir(File resultRoot, String sampleBaseName) {
		return new File(searchResultsRoot(resultRoot), sampleBaseName);
	}

	/** Round-one (marker / genome discovery) folder for a sample. */
	public static File markerDir(File resultRoot, String sampleBaseName) {
		return new File(sampleDir(resultRoot, sampleBaseName), MARKER_SEARCH_DIR);
	}

	/** Round-two (refined database) folder for a sample. */
	public static File sampleDbDir(File resultRoot, String sampleBaseName) {
		return new File(sampleDir(resultRoot, sampleBaseName), SAMPLE_DB_SEARCH_DIR);
	}

	/** Pooled validation output, alongside the other workflows' combined/. */
	public static File combinedDir(File resultRoot) {
		return new File(resultRoot, COMBINED_DIR);
	}

	/**
	 * Deliverable root: {@code combined/}, final tables, taxonomy, function,
	 * report, parameters and log. Same folder as {@code resultRoot} (the
	 * {@code opensearch_MAG_search} engine directory), matching FragPipe MAG.
	 * {@code engineMode} is ignored; closed+open MSFragger is the only mode.
	 */
	public static File modeRoot(File resultRoot, String engineMode) {
		migrateLegacyModeRoot(resultRoot);
		return resultRoot;
	}

	/**
	 * Legacy nested engine-mode folders. The cleaner still treats files inside
	 * them as final deliverables, and {@link MagResultLayoutUtil} does not
	 * migrate them into {@code search_results/} as if they were samples.
	 */
	public static List<String> modeFolderNames() {
		return Arrays.asList(LEGACY_MSFRAGGER_ONLY_DIR, LEGACY_MSFRAGGER_SAGE_DIR);
	}

	/**
	 * Lifts {@code msfragger_only/} up into the engine root when this folder is
	 * still on the pre-flatten layout. No-op when {@code combined/validated_psms.tsv}
	 * already exists here, or when there is no legacy folder. A leftover
	 * {@code msfragger_sage/} is not moved.
	 */
	public static void migrateLegacyModeRoot(File resultRoot) {
		if (resultRoot == null) {
			return;
		}
		File legacy = new File(resultRoot, LEGACY_MSFRAGGER_ONLY_DIR);
		if (!legacy.isDirectory()) {
			return;
		}
		if (isSatisfied(validatedPsms(resultRoot))) {
			return;
		}
		File[] children = legacy.listFiles();
		if (children == null) {
			return;
		}
		for (File child : children) {
			File dest = new File(resultRoot, child.getName());
			if (dest.exists()) {
				LOGGER.warn("Legacy lift skipped for {} - {} already exists", child.getAbsolutePath(),
						dest.getAbsolutePath());
				continue;
			}
			try {
				Files.move(child.toPath(), dest.toPath(), StandardCopyOption.ATOMIC_MOVE);
				LOGGER.info("Lifted open-search result {} -> {}", child.getAbsolutePath(), dest.getAbsolutePath());
			} catch (IOException atomicMoveFailed) {
				try {
					Files.move(child.toPath(), dest.toPath());
					LOGGER.info("Lifted open-search result {} -> {} (non-atomic fallback)", child.getAbsolutePath(),
							dest.getAbsolutePath());
				} catch (IOException e2) {
					LOGGER.error("Failed to lift {} -> {}: {}", child.getAbsolutePath(), dest.getAbsolutePath(),
							e2.getMessage(), e2);
				}
			}
		}
		String[] leftover = legacy.list();
		if (leftover != null && leftover.length == 0 && !legacy.delete()) {
			LOGGER.warn("Could not remove empty leftover {}", legacy.getAbsolutePath());
		}
	}

	/**
	 * Pooled validation output: {@code combined/} inside {@link #modeRoot}.
	 */
	public static File combinedDir(File resultRoot, String engineMode) {
		return combinedDir(modeRoot(resultRoot, engineMode));
	}

	/**
	 * Filesystem-safe key for an engine mode, e.g. "MSFragger + Sage" becomes
	 * {@code msfragger_sage}. Kept deterministic so a rerun of the same mode
	 * resolves to the same folder and the resume contract still holds.
	 */
	public static String engineKey(String engineMode) {
		if (engineMode == null || engineMode.isBlank()) {
			return "msfragger";
		}
		String key = engineMode.toLowerCase(java.util.Locale.ROOT).replaceAll("[^a-z0-9]+", "_")
				.replaceAll("^_+|_+$", "");
		return key.isEmpty() ? "msfragger" : key;
	}

	/**
	 * Round-two search folder for a sample, scoped by engine.
	 * <p>
	 * MSFragger keeps the plain {@code sample_db/} it has always used, so
	 * existing result folders stay readable and its searches are reused
	 * unchanged when a second engine is added. Any other engine gets its own
	 * suffixed folder rather than colliding on the one pepXML path per stage.
	 */
	public static File sampleDbDir(File resultRoot, String sampleBaseName, String engineName) {
		if (engineName == null || engineName.isBlank() || "msfragger".equals(engineKey(engineName))) {
			return sampleDbDir(resultRoot, sampleBaseName);
		}
		return new File(sampleDir(resultRoot, sampleBaseName),
				SAMPLE_DB_SEARCH_DIR + "_" + engineKey(engineName));
	}

	// ---------------------------------------------------------- sentinels

	/** Shared target/decoy marker database for round one. */
	/**
	 * Fallback only. The marker target-decoy FASTA normally lives beside
	 * {@code Marker4.fasta} in the MAG database folder, resolved through
	 * {@code MagDbItem.getHapFastaTDFile(selectorMode)} — the same path the
	 * FragPipe MAG workflow uses, so both search one shared file and neither
	 * rebuilds MSFragger's fragment index unnecessarily. This result-folder path
	 * is used only when the MAG database cannot supply one.
	 */
	public static File markerTargetDecoyFasta(File resultRoot) {
		return new File(resultRoot, "marker_target_decoy.fasta");
	}

	/**
	 * Round-one search sentinel. MSFragger writes pepXML next to its input, so
	 * it is copied into the sample folder once the search completes — which also
	 * makes it the last file the step writes, satisfying rule 1.
	 */
	public static File markerPepXml(File resultRoot, String sampleBaseName) {
		return new File(markerDir(resultRoot, sampleBaseName), sampleBaseName + ".pepXML");
	}

	/**
	 * Round-one MSFragger parameters, beside {@link #markerPepXml}.
	 * <p>
	 * Round one is one batch invocation over every sample, so a single parameter
	 * file governs them all; each sample gets its own copy rather than a shared
	 * one in a folder of its own, so that every search folder stands alone.
	 */
	public static File markerParams(File resultRoot, String sampleBaseName) {
		return new File(markerDir(resultRoot, sampleBaseName), SEARCH_PARAMS_FILE);
	}

	/** Round-two MSFragger parameters, beside {@link #sampleDbPepXml}. */
	public static File sampleDbParams(File resultRoot, String sampleBaseName) {
		return new File(sampleDbDir(resultRoot, sampleBaseName), SEARCH_PARAMS_FILE);
	}

	/**
	 * Round-one unmodified PSM set, consumed by genome selection. Written after
	 * the pepXML, so it gates zero-peak extraction separately from the search.
	 */
	public static File unmodifiedPsms(File resultRoot, String sampleBaseName) {
		return new File(markerDir(resultRoot, sampleBaseName), "unmodified_psms.tsv");
	}

	/**
	 * Written when stage-1 evidence is too thin to build a refined database for
	 * this sample. The sample then skips stage-2 identification but can still
	 * enter FlashLFQ quantification (match-between-runs) with the rest of the
	 * cohort, matching DIA Deep per-file's continue-with-successful gate.
	 */
	public static File identificationSkipped(File resultRoot, String sampleBaseName) {
		return new File(markerDir(resultRoot, sampleBaseName), "identification_skipped.txt");
	}

	/**
	 * The {@code .pin} MSFragger writes beside the marker pepXML, and the file
	 * handed to MSBooster unchanged.
	 * <p>
	 * It is named {@code <sample>.pin} and must stay that way: MSBooster pairs a
	 * pin with its spectra file by base name ({@code mainsteps.PinMzmlMatcher}),
	 * so any suffix sends it looking for {@code <sample>_suffix.mzML} and it
	 * aborts with "No matching mzml/mgf file" - an error that names the mzML
	 * directory and so points at the wrong thing entirely. This mirrors FragPipe,
	 * which feeds MSBooster its own {@code <run>/<run>.pin} as it stands.
	 */
	public static File markerPin(File resultRoot, String sampleBaseName) {
		return new File(markerDir(resultRoot, sampleBaseName), sampleBaseName + ".pin");
	}

	/**
	 * What MSBooster writes beside {@link #markerPin} under
	 * {@code renamePin = 1}. The name is MSBooster's to choose, so it is derived
	 * from its own definition rather than spelled out again here.
	 */
	public static File markerEditedPin(File resultRoot, String sampleBaseName) {
		return MSBoosterHeadlessRunner.editedPin(markerPin(resultRoot, sampleBaseName));
	}

	/**
	 * MSBooster's parameter file for the stage-1 batch, at the search-results
	 * root rather than in any one sample's folder.
	 * <p>
	 * One invocation covers every sample - the marker search is itself a single
	 * batched MSFragger run - so the parameters govern the batch, not a sample,
	 * and putting them in the first sample's folder would misdescribe them. This
	 * is the same shape FragPipe uses: {@code msbooster_params.txt} beside the
	 * per-run folders it names, not inside one of them.
	 */
	public static File msboosterParams(File resultRoot) {
		return new File(searchResultsRoot(resultRoot), MSBoosterHeadlessRunner.PARAMS_FILE);
	}

	/**
	 * Where a stage's MSBooster predictions are kept: {@code MSBooster_marker} or
	 * {@code MSBooster_sample_db} at the search-results root.
	 * <p>
	 * Stage-named, because MSBooster is not: it always writes to a folder called
	 * {@code MSBooster} under the deepest common directory of the pins it was
	 * given, and a stage-one batch over {@code <sample>/marker/*.pin} and a
	 * stage-two batch over {@code <sample>/sample_db/*.pin} both resolve to the
	 * same parent. Left alone, the second run would overwrite the first's
	 * predictions with no warning.
	 *
	 * @param searchDir {@link #MARKER_SEARCH_DIR} or
	 *                  {@link #SAMPLE_DB_SEARCH_DIR}.
	 */
	public static File msboosterDir(File resultRoot, String searchDir) {
		return new File(searchResultsRoot(resultRoot), "MSBooster_" + searchDir);
	}


	/**
	 * Stage-1 expected-shift audit (open marker search only): every UNIMOD shift with
	 * PSMs in its window for this sample, admitted or not, with the numbers behind each decision. Beside the
	 * PSM file it explains, so a marker folder stays self-describing.
	 */
	public static File expectedPeaks(File resultRoot, String sampleBaseName) {
		return new File(markerDir(resultRoot, sampleBaseName), "expected_peaks.tsv");
	}

	/**
	 * Refined, genome-specific database for this sample (targets + 1:1 decoys).
	 * Named {@code <sample>.fasta} in the marker folder, the same place and
	 * shape the other MAG engines write their sample-specific database.
	 */
	public static File refinedFasta(File resultRoot, String sampleBaseName) {
		return new File(markerDir(resultRoot, sampleBaseName), sampleBaseName + ".fasta");
	}

	/** Round-two search sentinel. */
	public static File sampleDbPepXml(File resultRoot, String sampleBaseName) {
		return new File(sampleDbDir(resultRoot, sampleBaseName), sampleBaseName + ".pepXML");
	}

	/**
	 * Closed search of the refined database, the first stage of the closed-then-open
	 * identification ({@code DDAMagOpenSearchParameters#isClosedStageEnabled()}). Its own folder
	 * beside {@code sample_db/}, so the open search and its Crystal-C output are untouched.
	 */
	public static File sampleDbClosedDir(File resultRoot, String sampleBaseName) {
		return new File(sampleDir(resultRoot, sampleBaseName), SAMPLE_DB_CLOSED_SEARCH_DIR);
	}

	/** Closed-stage search sentinel. */
	public static File sampleDbClosedPepXml(File resultRoot, String sampleBaseName) {
		return new File(sampleDbClosedDir(resultRoot, sampleBaseName), sampleBaseName + ".pepXML");
	}

	/** Closed-stage MSFragger parameters, beside {@link #sampleDbClosedPepXml}. */
	public static File sampleDbClosedParams(File resultRoot, String sampleBaseName) {
		return new File(sampleDbClosedDir(resultRoot, sampleBaseName), SEARCH_PARAMS_FILE);
	}

	/**
	 * Closed-then-open run summary: stage and class counts, peptide-level thresholds, the
	 * closed/open disagreement and the per-substitution-peak closed agreement.
	 */
	public static File closedOpenSummary(File resultRoot, String engineMode) {
		return new File(combinedDir(resultRoot, engineMode), CLOSED_OPEN_SUMMARY_FILE);
	}

	/** Joint closed+open identification: group counts, peptide walks, model coefficients. */
	public static File jointSummary(File resultRoot, String engineMode) {
		return new File(combinedDir(resultRoot, engineMode), "joint_identification_summary.tsv");
	}

	/** Joint identification: the per-peak fragment test ({@link PeakFragmentTest}), one row per established peak. */
	public static File peakFragmentTest(File resultRoot, String engineMode) {
		return new File(combinedDir(resultRoot, engineMode), "peak_fragment_test.tsv");
	}

	/** Joint identification: target/decoy AUC of every feature, per FDR group. */
	public static File jointFeatureReport(File resultRoot, String engineMode) {
		return new File(combinedDir(resultRoot, engineMode), "joint_feature_report.tsv");
	}

	/** Joint identification: DIA-NN prediction inputs and outputs, per run. */
	public static File jointPredictionDir(File resultRoot, String engineMode) {
		return new File(combinedDir(resultRoot, engineMode), "joint_predictions");
	}

	/**
	 * Closed-only, open-only, cascade and joint identifications of the same spectra, one row per
	 * spectrum winner with q &lt;= 0.05 per method, for comparing the methods at matched error.
	 */
	public static File identificationArms(File resultRoot, String engineMode) {
		return new File(combinedDir(resultRoot, engineMode), "identification_arms.tsv.gz");
	}

	/** Proteins of every peptide in {@link #identificationArms}. */
	public static File identificationArmPeptides(File resultRoot, String engineMode) {
		return new File(combinedDir(resultRoot, engineMode), "identification_arm_peptides.tsv.gz");
	}

	/** Decoy-estimated counts per method at 1%. */
	public static File identificationArmsSummary(File resultRoot, String engineMode) {
		return new File(combinedDir(resultRoot, engineMode), "identification_arms_summary.tsv");
	}

	/**
	 * Which identification scheme the validation outputs beside it were produced under (open only,
	 * or closed then open). Rewritten only when the scheme changes, so it is an input of pooled
	 * validation that is newer than its outputs exactly when the switch was flipped.
	 */
	public static File identificationMode(File resultRoot, String engineMode) {
		return new File(combinedDir(resultRoot, engineMode), IDENTIFICATION_MODE_FILE);
	}

	public static final String IDENTIFICATION_MODE_FILE = "identification_mode.txt";
	public static final String CLOSED_OPEN_SUMMARY_FILE = "closed_open_summary.tsv";

	/**
	 * Round one after Crystal-C: the precursor-refined marker pepXML, and the
	 * sentinel for the stage-1 refinement step.
	 * <p>
	 * Stage 1 admits only the zero family plus tested UNIMOD shifts, so
	 * it is the step the refinement helps most: a PSM whose precursor was
	 * mis-picked sits a few Da off zero and is discarded, however good the match
	 * is. Measured on 12mix (2026-09-14) at the stage-2 stage, where the same
	 * correction was run first: rank-1 target PSMs inside |delta| &lt; 0.05 Da went
	 * from 310,049 to 581,451 at 1% FDR, at unchanged purity against the 12
	 * spiked genomes (96.98% -> 97.03%).
	 */
	public static File markerCorrectedPepXml(File resultRoot, String sampleBaseName) {
		return CrystalCHeadlessRunner.correctedPepXml(markerDir(resultRoot, sampleBaseName),
				markerPepXml(resultRoot, sampleBaseName));
	}

	/**
	 * Round two after Crystal-C: the precursor-refined pepXML, and the sentinel
	 * for the refinement step.
	 * <p>
	 * The name is Crystal-C's to choose, so it is derived from its own definition
	 * rather than spelled out again here. Everything downstream - the pooled
	 * reader, delta-mass fitting, the classifier - consumes this file rather than
	 * {@link #sampleDbPepXml} when the step is enabled; the uncorrected pepXML is
	 * kept beside it as the input this was derived from, and remains the search
	 * step's own sentinel.
	 */
	public static File sampleDbCorrectedPepXml(File resultRoot, String sampleBaseName) {
		return CrystalCHeadlessRunner.correctedPepXml(sampleDbDir(resultRoot, sampleBaseName),
				sampleDbPepXml(resultRoot, sampleBaseName));
	}

	/**
	 * Crystal-C's per-PSM report: one row per corrected spectrum with the old and
	 * new precursor mass, charge and delta mass, and its reason for each change.
	 * A diagnostic artifact worth keeping - it is the only place the refinement
	 * explains itself.
	 */
	public static File crystalcReport(File resultRoot, String sampleBaseName) {
		return CrystalCHeadlessRunner.report(sampleDbDir(resultRoot, sampleBaseName),
				sampleDbPepXml(resultRoot, sampleBaseName));
	}

	/** Crystal-C's parameter file for this sample, kept for provenance. */
	public static File crystalcParams(File resultRoot, String sampleBaseName) {
		return new File(sampleDbDir(resultRoot, sampleBaseName),
				"crystalc-" + sampleDbPepXml(resultRoot, sampleBaseName).getName() + ".params");
	}

	/** Where each screened peak's localised PSMs sit; see {@code PeakSiteShare}. */
	public static File peakSiteShares(File resultRoot, String engineMode) {
		return new File(combinedDir(resultRoot, engineMode), "peak_site_shares.tsv");
	}

	/**
	 * Per-peak audit of the substitution analysis: which fitted peaks sit at a
	 * single amino-acid substitution mass, and how many of their PSMs the
	 * localization evidence actually agrees with.
	 */
	public static File substitutionPeaks(File resultRoot, String engineMode) {
		return inSequenceVariants(resultRoot, engineMode, "substitution_peaks.tsv");
	}

	/**
	 * Proposed sequences for the bounded re-search: each assigned peptide with
	 * the one residue changed that its delta mass and localization both support.
	 * <p>
	 * This is the only route to the majority case. Measured on 12mix, 70.5% of
	 * the PSMs on the +/-14.0157 Da peaks propose a sequence that is <em>absent
	 * from every database this run searched</em> - genome selection picked one
	 * strain and the peptide came from a relative - so no lookup can confirm
	 * them and only a search against these generated sequences can.
	 */
	public static File substitutionCandidates(File resultRoot, String engineMode) {
		return inSequenceVariants(resultRoot, engineMode, "substitution_candidates.fasta");
	}

	/**
	 * Site-level single-residue sequence variants called by
	 * {@link uk.ac.quadram.metapilot.core.opensearch.SequenceVariantCaller}: one
	 * row per protein site, with its evidence class, recurrence and whether the
	 * reference peptide was also observed. PSM membership is in
	 * {@link #substitutionForms(File, String)}'s {@code variant_id} column.
	 */
	public static File sequenceVariants(File resultRoot, String engineMode) {
		return inSequenceVariants(resultRoot, engineMode, "sequence_variants.tsv");
	}

	/**
	 * Reference proteins with the recurrent variants applied, one entry per
	 * (variant, protein), for a follow-up search that should identify the variant
	 * peptides as ordinary unmodified peptides.
	 */
	public static File sequenceVariantsFasta(File resultRoot, String engineMode) {
		return inSequenceVariants(resultRoot, engineMode, "sequence_variants.fasta");
	}

	/** {@code sequence_variants/} under the engine result root. */
	public static File sequenceVariantsDir(File resultRoot, String engineMode) {
		return new File(modeRoot(resultRoot, engineMode), SEQUENCE_VARIANTS_DIR);
	}

	/**
	 * FlashLFQ identification input for substitution PSMs that are not in
	 * {@code validated_psms.tsv}, plus the validated unmodified PSMs of the same base
	 * peptides as references ({@code quant_role}). Match-between-runs and
	 * normalisation are off for the run that reads it, and nothing outside
	 * {@code sequence_variants/} reads the result.
	 */
	public static final String SUBSTITUTION_QUANT_PSMS_FILE = "substitution_psms.tsv";

	/**
	 * FlashLFQ's peptide table from {@link #substitutionQuantPsms}, under its own name,
	 * with role, variant ids, variant/reference ratios and shared main-run features appended.
	 */
	public static final String SUBSTITUTION_QUANTIFIED_PEPTIDES_FILE = "substitution_quantified_peptides.tsv";

	/** FlashLFQ's peak table from the same run, under its own name. */
	public static final String SUBSTITUTION_QUANTIFIED_PEAKS_FILE = "substitution_quantified_peaks.tsv";

	public static File substitutionQuantPsms(File resultRoot, String engineMode) {
		return inSequenceVariants(resultRoot, engineMode, SUBSTITUTION_QUANT_PSMS_FILE);
	}

	public static File substitutionQuantifiedPeptides(File resultRoot, String engineMode) {
		return inSequenceVariants(resultRoot, engineMode, SUBSTITUTION_QUANTIFIED_PEPTIDES_FILE);
	}

	public static File substitutionQuantifiedPeaks(File resultRoot, String engineMode) {
		return inSequenceVariants(resultRoot, engineMode, SUBSTITUTION_QUANTIFIED_PEAKS_FILE);
	}

	/** FlashLFQ's peak table from quantifying {@code validated_psms.tsv}. */
	public static File quantifiedPeaks(File resultRoot, String engineMode) {
		return new File(combinedDir(resultRoot, engineMode), "QuantifiedPeaks.tsv");
	}

	/**
	 * A file in {@link #SEQUENCE_VARIANTS_DIR}. A copy left in {@code combined/}
	 * by an earlier run is moved here when this path is still empty, so a
	 * finished result is not treated as missing.
	 */
	private static File inSequenceVariants(File resultRoot, String engineMode, String name) {
		File current = new File(sequenceVariantsDir(resultRoot, engineMode), name);
		File legacy = new File(combinedDir(resultRoot, engineMode), name);
		if (current.isFile() || !legacy.isFile()) {
			return current;
		}
		File dir = current.getParentFile();
		if (dir != null && (dir.isDirectory() || dir.mkdirs())) {
			try {
				Files.move(legacy.toPath(), current.toPath(), StandardCopyOption.ATOMIC_MOVE);
				LOGGER.info("Moved open-search variant output {} -> {}", legacy.getAbsolutePath(),
						current.getAbsolutePath());
				return current;
			} catch (IOException atomicMoveFailed) {
				try {
					Files.move(legacy.toPath(), current.toPath());
					LOGGER.info("Moved open-search variant output {} -> {} (non-atomic fallback)",
							legacy.getAbsolutePath(), current.getAbsolutePath());
					return current;
				} catch (IOException e2) {
					LOGGER.warn("Could not move {} into {}: {}", legacy.getAbsolutePath(),
							current.getAbsolutePath(), e2.getMessage());
				}
			}
		}
		return legacy;
	}

	/** Fitted delta-mass peaks; a diagnostic artifact worth keeping. */
	public static File fittedPeaks(File resultRoot) {
		return new File(combinedDir(resultRoot), "delta_mass_peaks.tsv");
	}

	/** Engine-scoped; see {@link #combinedDir(File, String)}. */
	public static File fittedPeaks(File resultRoot, String engineMode) {
		return new File(combinedDir(resultRoot, engineMode), "delta_mass_peaks.tsv");
	}

	/** Classifier-scored PSMs, pre-validation. */
	public static File scoredPsms(File resultRoot) {
		return new File(combinedDir(resultRoot), "scored_psms.tsv");
	}

	/** Engine-scoped; see {@link #combinedDir(File, String)}. */
	public static File scoredPsms(File resultRoot, String engineMode) {
		return new File(combinedDir(resultRoot, engineMode), "scored_psms.tsv");
	}

	/** Validated PSMs. Sentinel for the whole validation stage. */
	public static File validatedPsms(File resultRoot) {
		return new File(combinedDir(resultRoot), "validated_psms.tsv");
	}

	/** Engine-scoped; see {@link #combinedDir(File, String)}. */
	public static File validatedPsms(File resultRoot, String engineMode) {
		return new File(combinedDir(resultRoot, engineMode), "validated_psms.tsv");
	}

	/** PSMs supported by an accepted fitted delta-mass hypothesis. */
	public static File peakValidatedPsms(File resultRoot) {
		return new File(combinedDir(resultRoot), "peak_validated_psms.tsv");
	}

	/** Engine-scoped; see {@link #combinedDir(File, String)}. */
	public static File peakValidatedPsms(File resultRoot, String engineMode) {
		return new File(combinedDir(resultRoot, engineMode), "peak_validated_psms.tsv");
	}

	/** {@code substitution_forms.tsv}: one row per flagged PSM and its candidate form. */
	public static final String SUBSTITUTION_FORMS_FILE = "substitution_forms.tsv";

	/**
	 * The substitution-candidate relations for {@code validated_psms.tsv}: one row
	 * per PSM that {@link uk.ac.quadram.metapilot.core.opensearch.SubstitutionReassigner}
	 * flagged, carrying every candidate form and the verdict. Kept out of
	 * {@code validated_psms.tsv} itself so that file holds exactly one row per PSM
	 * and can be handed to FlashLFQ unmodified.
	 * <p>
	 * New runs write this under {@code sequence_variants/}. A copy already in
	 * {@code combined/} is used as it stands and is not moved.
	 */
	public static File substitutionForms(File resultRoot) {
		return substitutionForms(resultRoot, null);
	}

	/** See {@link #substitutionForms(File)}. */
	public static File substitutionForms(File resultRoot, String engineMode) {
		File current = substitutionFormsOutput(resultRoot, engineMode);
		if (current.isFile() && current.length() > 0) {
			return current;
		}
		File legacy = new File(combinedDir(resultRoot, engineMode), SUBSTITUTION_FORMS_FILE);
		if (legacy.isFile() && legacy.length() > 0) {
			return legacy;
		}
		return current;
	}

	/** Where a new run writes {@code substitution_forms.tsv}. */
	public static File substitutionFormsOutput(File resultRoot, String engineMode) {
		return new File(sequenceVariantsDir(resultRoot, engineMode), SUBSTITUTION_FORMS_FILE);
	}

	/**
	 * FlashLFQ's own output from quantifying {@link #validatedPsms(File)}: a
	 * peptide-level table with one row per distinct peptide and its "Protein
	 * Groups" resolved from every PSM for that peptide across the single
	 * combined run over every sample. The filename is FlashLFQ's own fixed
	 * output name, not one this workflow chooses - it always lands beside the
	 * {@code --idt} file it was given, so this is
	 * {@link #combinedDir(File)}/QuantifiedPeptides.tsv, matching the shared
	 * {@code combined/QuantifiedPeptides.tsv} convention every other DDA MAG
	 * engine (pFind, FragPipe, Sage) already uses - see
	 * {@code MagResultLayoutUtil}'s layout diagram.
	 */
	public static File quantifiedPeptides(File resultRoot) {
		return new File(combinedDir(resultRoot), "QuantifiedPeptides.tsv");
	}

	/** Engine-scoped; see {@link #combinedDir(File, String)}. */
	public static File quantifiedPeptides(File resultRoot, String engineMode) {
		return new File(combinedDir(resultRoot, engineMode), "QuantifiedPeptides.tsv");
	}

	/** Confident peptide IDs whose observed mass is not assigned to a fitted peak. */
	public static File unexplainedPsms(File resultRoot) {
		return new File(combinedDir(resultRoot), "unexplained_psms.tsv");
	}

	/** Engine-scoped; see {@link #combinedDir(File, String)}. */
	public static File unexplainedPsms(File resultRoot, String engineMode) {
		return new File(combinedDir(resultRoot, engineMode), "unexplained_psms.tsv");
	}

	/** Per-peak gate verdicts: counts, raw and BH-adjusted p-values, retained PSMs. */
	public static File peakVerdicts(File resultRoot) {
		return new File(combinedDir(resultRoot), "modification_verdicts.tsv");
	}

	/** Engine-scoped; see {@link #combinedDir(File, String)}. */
	public static File peakVerdicts(File resultRoot, String engineMode) {
		return new File(combinedDir(resultRoot, engineMode), "modification_verdicts.tsv");
	}

	/**
	 * One row per <em>accepted</em> delta-mass peak, joining the fitted-peak
	 * statistics and FDR verdict with candidate explanations for its mass: every
	 * bundled chemical modification and single amino-acid substitution whose mass
	 * falls within tolerance of the peak's centroid. The "what might this delta
	 * mass be" companion to {@link #fittedPeaks} and {@link #peakVerdicts}, which
	 * say whether a peak is real but not what it means.
	 */
	public static File potentialModifications(File resultRoot) {
		return new File(combinedDir(resultRoot), "potential_modifications.tsv");
	}

	/** Engine-scoped; see {@link #combinedDir(File, String)}. */
	public static File potentialModifications(File resultRoot, String engineMode) {
		return new File(combinedDir(resultRoot, engineMode), "potential_modifications.tsv");
	}

	/**
	 * Pangenome variant check, one row per substitution candidate (PSM x site)
	 * with its UHGP-100 lookup and class. Written only when the add-on runs;
	 * see {@code uk.ac.quadram.metapilot.core.variant.PangenomeVariantCheck}.
	 */
	public static File variantCandidates(File resultRoot, String engineMode) {
		return inSequenceVariants(resultRoot, engineMode,
				uk.ac.quadram.metapilot.core.variant.PangenomeVariantCheck.CANDIDATES_FILE);
	}

	/** {@link #variantCandidates} rolled up to one row per protein site and substitution. */
	public static File variantSites(File resultRoot, String engineMode) {
		return inSequenceVariants(resultRoot, engineMode,
				uk.ac.quadram.metapilot.core.variant.PangenomeVariantCheck.SITES_FILE);
	}

	/** Per-substitution / per-peak summary of {@link #variantCandidates}, against matched random-substitution controls. */
	public static File variantPangenomeSummary(File resultRoot, String engineMode) {
		return inSequenceVariants(resultRoot, engineMode,
				uk.ac.quadram.metapilot.core.variant.PangenomeVariantCheck.SUMMARY_FILE);
	}

	/** Per-species rollup of {@link #variantCandidates}, written last by the add-on. */
	public static File variantSpeciesSummary(File resultRoot, String engineMode) {
		return inSequenceVariants(resultRoot, engineMode,
				uk.ac.quadram.metapilot.core.variant.PangenomeVariantCheck.SPECIES_FILE);
	}

	/** Column dictionary written beside the other UHGP tables. */
	public static File variantColumns(File resultRoot, String engineMode) {
		return inSequenceVariants(resultRoot, engineMode,
				uk.ac.quadram.metapilot.core.variant.PangenomeVariantCheck.COLUMNS_FILE);
	}

	// --- deliverables -----------------------------------------------------
	// Same filenames the other MAG workflows produce, so downstream report
	// generation and the cleaner's ROOT_FINAL_FILES rule pick them up without
	// workflow-specific handling. They live in a {@link #modeRoot}, which is
	// what these take.

	/** Final peptide table in a {@link #modeRoot}. */
	public static File finalPeptides(File modeRoot) {
		return new File(modeRoot, "final_peptides.tsv");
	}

	/** Razor protein table in a {@link #modeRoot}. */
	public static File razorProteins(File modeRoot) {
		return new File(modeRoot, "razor_proteins.tsv");
	}

	/** All-protein table in a {@link #modeRoot}. */
	public static File allProteins(File modeRoot) {
		return new File(modeRoot, "all_proteins.tsv");
	}

	/**
	 * Genome-selection audit artifacts, mirroring the FragPipe workflow's
	 * {@code *_evidence/_ranked/_selected.tsv} convention so the cleaner's
	 * existing keep-rule shape carries over.
	 */
	public static List<String> genomeSelectionAuditSuffixes() {
		return Arrays.asList("_evidence.tsv", "_ranked.tsv", "_selected.tsv");
	}

	/**
	 * Files that must survive even a deep clean, because a later run's resume
	 * check would otherwise skip the step that produces them and never
	 * regenerate them. Paths are relative to the result root.
	 */
	public static List<String> deepCleanKeepFileNames() {
		// SEARCH_PARAMS_FILE is kept for the reason in rule 2 above, applied to
		// provenance rather than data: the pepXML beside it is a resume sentinel,
		// so deleting the parameters leaves a search that will be skipped on the
		// next run with no record of what produced it. That is precisely the
		// silent drift this class exists to prevent.
		return Arrays.asList("unmodified_psms.tsv", "scored_psms.tsv", "validated_psms.tsv", "peak_validated_psms.tsv",
				"unexplained_psms.tsv", "modification_verdicts.tsv", "delta_mass_peaks.tsv",
				"marker_target_decoy.fasta", "final_peptides.tsv", "razor_proteins.tsv", "all_proteins.tsv",
				SEARCH_PARAMS_FILE, IDENTIFICATION_MODE_FILE, CLOSED_OPEN_SUMMARY_FILE,
				// Regenerating these means re-running validation (the reassigner's
				// result is not on disk), so they are kept like its sentinels.
				uk.ac.quadram.metapilot.core.variant.PangenomeVariantCheck.CANDIDATES_FILE,
				uk.ac.quadram.metapilot.core.variant.PangenomeVariantCheck.SITES_FILE,
				uk.ac.quadram.metapilot.core.variant.PangenomeVariantCheck.SUMMARY_FILE,
				uk.ac.quadram.metapilot.core.variant.PangenomeVariantCheck.SPECIES_FILE,
				uk.ac.quadram.metapilot.core.variant.PangenomeVariantCheck.COLUMNS_FILE,
				SUBSTITUTION_QUANT_PSMS_FILE, SUBSTITUTION_QUANTIFIED_PEPTIDES_FILE, SUBSTITUTION_QUANTIFIED_PEAKS_FILE,
				SUBSTITUTION_FORMS_FILE);
	}

	/** True when a sentinel counts as satisfied. */
	public static boolean isSatisfied(File sentinel) {
		return sentinel != null && sentinel.exists() && sentinel.length() > 0;
	}

	/**
	 * True when a sentinel is satisfied <em>and</em> no input it was derived from
	 * is newer than it.
	 * <p>
	 * Existence alone is not enough for a step whose inputs can change without
	 * the step being asked to re-run. Switching stage-1 rescoring on is exactly
	 * that: MSBooster writes a new edited pin, but the refined FASTA from the
	 * previous run still exists, so genome selection is skipped and stage 2
	 * silently searches a database built by the cut the rescoring was meant to
	 * replace. Nothing in the output says so. This is the make-style check that
	 * catches it - an input newer than the artefact means the artefact is stale.
	 * <p>
	 * Inputs that do not exist are ignored rather than treated as changes, so an
	 * optional step that never ran does not force everything downstream to
	 * re-run.
	 *
	 * @return false when the sentinel is missing, empty, or older than any
	 *         existing input.
	 */
	public static boolean isSatisfied(File sentinel, File... inputs) {
		if (!isSatisfied(sentinel)) {
			return false;
		}
		if (inputs == null) {
			return true;
		}
		for (File input : inputs) {
			if (input != null && input.exists() && input.lastModified() > sentinel.lastModified()) {
				return false;
			}
		}
		return true;
	}

	/** The newest existing input among {@code inputs}, or null when none exist. */
	public static File newestInput(File... inputs) {
		File newest = null;
		if (inputs == null) {
			return null;
		}
		for (File input : inputs) {
			if (input != null && input.exists() && (newest == null || input.lastModified() > newest.lastModified())) {
				newest = input;
			}
		}
		return newest;
	}
}
