package uk.ac.quadram.metapilot.core.runner;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.nio.file.Files;
import java.nio.file.StandardCopyOption;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.EnumMap;
import java.util.HashMap;
import java.util.HashSet;
import java.util.IdentityHashMap;
import java.util.TreeMap;
import java.util.TreeSet;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.concurrent.atomic.AtomicReference;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import uk.ac.quadram.metapilot.core.crystalc.CrystalCHeadlessRunner;
import uk.ac.quadram.metapilot.core.mag.DdaMagProteinInferenceService;
import uk.ac.quadram.metapilot.core.mag.GenomeEvidenceTest;
import uk.ac.quadram.metapilot.core.mag.MagGenomeSelectionService;
import uk.ac.quadram.metapilot.core.mag.MagResultLayoutUtil;
import uk.ac.quadram.metapilot.core.msbooster.MSBoosterHeadlessRunner;
import uk.ac.quadram.metapilot.core.msbooster.MSBoosterPin;
import uk.ac.quadram.metapilot.core.opensearch.Stage1Rescorer;
import uk.ac.quadram.metapilot.core.msfragger.MSFraggerHeadlessRunner;
import uk.ac.quadram.metapilot.core.opensearch.ClosedOpenCascade;
import uk.ac.quadram.metapilot.core.opensearch.IdentificationArms;
import uk.ac.quadram.metapilot.core.opensearch.ShiftLocaliser;
import uk.ac.quadram.metapilot.core.opensearch.XGBoostTreeClassifier;
import uk.ac.quadram.metapilot.core.opensearch.JointCandidate;
import uk.ac.quadram.metapilot.core.opensearch.JointIdentification;
import uk.ac.quadram.metapilot.core.opensearch.PredictedSpectrumFeatures;
import uk.ac.quadram.metapilot.core.opensearch.DeltaMassEvidenceGate;
import uk.ac.quadram.metapilot.core.opensearch.SequenceVariantCaller;
import uk.ac.quadram.metapilot.core.opensearch.SubstitutionReassigner;
import uk.ac.quadram.metapilot.core.opensearch.TrypticPeptideIndex;
import uk.ac.quadram.metapilot.core.opensearch.DeltaMassPeakFitter;
import uk.ac.quadram.metapilot.core.opensearch.ExpectedPeakSelector;
import uk.ac.quadram.metapilot.core.opensearch.ModificationMassDatabase;
import uk.ac.quadram.metapilot.core.opensearch.ModificationValidator;
import uk.ac.quadram.metapilot.core.opensearch.OpenSearchFeatures;
import uk.ac.quadram.metapilot.core.opensearch.OpenSearchLayout;
import uk.ac.quadram.metapilot.core.opensearch.PeakFragmentTest;
import uk.ac.quadram.metapilot.core.opensearch.PeakSiteShare;
import uk.ac.quadram.metapilot.core.opensearch.OpenSearchSpectra;
import uk.ac.quadram.metapilot.core.opensearch.OpenSearchEntrapment;
import uk.ac.quadram.metapilot.core.opensearch.OpenSearchMetaPsm;
import uk.ac.quadram.metapilot.core.opensearch.OpenSearchPepxmlReader;
import uk.ac.quadram.metapilot.core.opensearch.OpenSearchPsm;
import uk.ac.quadram.metapilot.core.opensearch.OpenSearchValidatedPsmReader;
import uk.ac.quadram.metapilot.core.opensearch.SemiSupervisedPsmClassifier;
import uk.ac.quadram.metapilot.core.opensearch.SpectrumLevelCollapse;
import uk.ac.quadram.metapilot.core.opensearch.SubstitutionMassTable;
import uk.ac.quadram.metapilot.core.opensearch.TargetDecoyFdr;
import uk.ac.quadram.metapilot.core.opensearch.UnmodifiedPsmSelector;
import uk.ac.quadram.metapilot.core.opensearch.CombinedPsmCandidate;
import uk.ac.quadram.metapilot.core.opensearch.CombinedPsmClassifier;
import uk.ac.quadram.metapilot.core.opensearch.PerEngineClassifier;
import uk.ac.quadram.metapilot.core.prodb.TargetDecoyFasta;
import uk.ac.quadram.metapilot.core.sage.SageOpenSearchReader;
import uk.ac.quadram.metapilot.core.sage.SageOpenSearchRunner;
import uk.ac.quadram.metapilot.core.variant.OpenSearchVariantCandidateSource;
import uk.ac.quadram.metapilot.core.variant.PangenomeVariantCheck;
import uk.ac.quadram.metapilot.core.variant.VariantCandidate;
import uk.ac.quadram.metapilot.core.variant.VariantFormRanking;
import uk.ac.quadram.metapilot.core.workflow.ProgressCallback;
import uk.ac.quadram.metapilot.quant.flashLFQ.FlashLfqTask;
import uk.ac.quadram.metapilot.task.dda.mag.par.DDAMagOpenSearchParameterIO;
import uk.ac.quadram.metapilot.task.dda.mag.par.DDAMagOpenSearchParameters;
import uk.ac.quadram.metapilot.task.io.pep.MetaPSM;
import uk.ac.quadram.metapilot.task.io.qpx.QpxPsmParquetWriter;
import uk.ac.quadram.metapilot.task.mag.GenomeSelectionCut;
import uk.ac.quadram.metapilot.task.mag.MagDbItem;
import uk.ac.quadram.metapilot.task.par.MetaConstants;
import uk.ac.quadram.metapilot.task.par.MetaData;

/**
 * Orchestrates the standalone MSFragger open-search MAG workflow.
 * <p>
 * Lives in {@code metapilot-core}, mirroring {@link MagFragPipeHeadlessRunner},
 * so the GUI, CLI and HPC entry points share one implementation rather than the
 * older pattern of burying logic in a GUI task.
 *
 * <h2>Steps</h2>
 * Step letters and messages follow {@link MagFragPipeHeadlessRunner}'s, and the
 * output folders are the shared DDA MAG ones, so a result folder from this
 * workflow reads the same as one from pFind, FragPipe or Sage.
 * <p>
 * Before Step A, Thermo {@code .raw} input is converted to mzML in the result
 * folder's {@code spectra/} and every step reads that; mzML input is read as
 * loaded.
 * See {@link RawSpectraConversion}.
 * <ol>
 * <li><b>A</b> — marker searches, per sample, against the marker/HAP database
 * (1:1 {@code rev_} decoys), written to
 * {@code search_results/<sample>/marker/}. Closed by default (ppm window,
 * isotope errors 0/1/2, Met oxidation and protein N-terminal acetylation; see
 * {@link #round1Settings}); {@code round1_search_mode = open} restores the
 * open search with Crystal-C. Either way only the zero family is used
 * downstream, plus, in open mode, PSMs at UNIMOD masses that pass
 * {@link ExpectedPeakSelector}'s criteria. All
 * samples share the database, so they go in one MSFragger invocation:
 * MSFragger reuses its fragment index across files in a single call, which a
 * per-sample loop would rebuild each time.</li>
 * <li><b>B</b> — select genomes from those unmodified PSMs and build the
 * sample-specific FASTA, exactly as the other MAG engines do. A sample with no
 * usable stage-1 targets can be skipped here (see
 * {@link PartialIdentificationHandler} / {@code --continue-with-successful-files})
 * instead of failing the whole cohort; that sample then has no stage-2 search
 * but still enters FlashLFQ with the others.</li>
 * <li><b>C</b> — sample-specific searches, per sample, into
 * {@code search_results/<sample>/sample_db/}. Batching is unavailable here
 * because each sample has its own database. Skipped stage-1 samples are omitted.</li>
 * <li><b>D</b> — pooled validation into {@code combined/}: delta-mass peak
 * fitting, one classifier over all PSMs, per-modification gate and local FDR.</li>
 * <li><b>F</b> — FlashLFQ quantification, a single combined run over the whole
 * {@code validated_psms.tsv} and every sample's raw/mzML directory together —
 * the same single-combined-run pattern {@code DDAMagPfindTask.quantify()}
 * uses, so FlashLFQ resolves cross-sample protein consistency the same way
 * for both workflows. Produces {@code combined/QuantifiedPeptides.tsv}.</li>
 * <li><b>G</b> — hand the validated output to the shared MAG post-analysis
 * services (genome scoring, protein inference, taxonomy, functional annotation,
 * reports), reused unchanged. The peptides {@link OpenSearchPostAnalysisTsvExporter}
 * builds are first rematched against step F's quantified output — see
 * {@code OpenSearchPostAnalysisTsvExporter#rematchWithFlashLfq} — so real
 * intensities and FlashLFQ's cross-sample protein sets reach genome scoring
 * and protein inference rather than phase-1's zero-intensity placeholders.</li>
 * </ol>
 *
 * <h2>Resume</h2>
 * Every expensive step is gated by a sentinel from {@link OpenSearchLayout},
 * which the cleaner reads from the same place so the two cannot drift.
 */
public final class MagOpenSearchHeadlessRunner {

	private static final Logger LOGGER = LogManager.getLogger(MagOpenSearchHeadlessRunner.class);

	// Progress bands, mirroring MagFragPipeHeadlessRunner's so the bar behaves
	// the same way across workflows.
	private static final double STEP_A_START_PERCENT = 5;
	private static final double STEP_A_END_PERCENT = 30;
	private static final double STEP_B_START_PERCENT = 30;
	private static final double STEP_B_END_PERCENT = 40;
	private static final double STEP_C_START_PERCENT = 40;
	private static final double STEP_C_END_PERCENT = 65;
	private static final double STEP_D_START_PERCENT = 65;
	private static final double STEP_D_END_PERCENT = 90;

	private MagOpenSearchHeadlessRunner() {
	}

	/**
	 * Optional GUI/interactive decision when some samples have no stage-1
	 * evidence for genome selection. Return {@code true} to continue stage-2
	 * identification + pooled validation using only the successful samples, then
	 * still quantify every raw file together in FlashLFQ; {@code false} to abort
	 * so the user can replace the failing files.
	 * <p>
	 * Mirrors {@code DiaMarkerDeepPerFileHeadlessRunner.PartialMarkerResultHandler}.
	 */
	@FunctionalInterface
	public interface PartialIdentificationHandler {
		/**
		 * @param skippedSampleNames display names of samples without a refined database
		 * @param successfulSampleCount number of samples that completed genome selection
		 * @return {@code true} to continue with successful samples; {@code false} to stop
		 */
		boolean continueWithSuccessfulFiles(List<String> skippedSampleNames, int successfulSampleCount)
				throws InterruptedException;
	}

	/** Opt-in for headless runs (CLI {@code --continue-with-successful-files}); default off. */
	public static final String CONTINUE_WITH_SUCCESSFUL_PROPERTY =
			"metapilot.mag_opensearch.continue_with_successful";

	static boolean continueWithSuccessfulRequested() {
		return Boolean.getBoolean(CONTINUE_WITH_SUCCESSFUL_PROPERTY);
	}

	public static final class Result {
		/** The engine search root, holding the shared {@code search_results/}. */
		public final File resultRoot;
		/** Same as {@link #resultRoot}: {@code opensearch_MAG_search/}. */
		public final File modeRoot;
		public final File validatedPsmsTsv;
		public final File peakVerdictsTsv;
		public final int validatedPsmCount;
		public final int acceptedModificationCount;

		Result(File resultRoot, File modeRoot, File validatedPsmsTsv, File peakVerdictsTsv, int validatedPsmCount,
				int acceptedModificationCount) {
			this.resultRoot = resultRoot;
			this.modeRoot = modeRoot;
			this.validatedPsmsTsv = validatedPsmsTsv;
			this.peakVerdictsTsv = peakVerdictsTsv;
			this.validatedPsmCount = validatedPsmCount;
			this.acceptedModificationCount = acceptedModificationCount;
		}
	}

	/**
	 * @param params        workflow parameters (the single source of truth for
	 *                      GUI, CLI and HPC alike)
	 * @param msfraggerJar  MSFragger.jar
	 * @param markerFasta   marker/HAP FASTA, targets only; decoys are added here
	 * @param cb            progress sink; may be null
	 * @deprecated no FlashLFQ executable is supplied, so Step F is skipped and
	 *             peptides keep phase-1 zero intensity — kept only so existing
	 *             callers keep compiling until they are updated to the overload
	 *             below. Prefer that overload for a real quantified run.
	 */
	@Deprecated
	public static Result run(DDAMagOpenSearchParameters params, File msfraggerJar, File markerFasta, MagDbItem mag,
			ProgressCallback cb) throws Exception {
		return run(params, msfraggerJar, markerFasta, mag, null, 1, cb);
	}

	/**
	 * @param params          workflow parameters (the single source of truth for
	 *                        GUI, CLI and HPC alike)
	 * @param msfraggerJar    MSFragger.jar
	 * @param markerFasta     marker/HAP FASTA, targets only; decoys are added here
	 * @param flashLfqExe     FlashLFQ CLI executable, sourced the same way callers
	 *                        already source {@code msfraggerJar} — this class does
	 *                        not know where tool installs live, only how to run
	 *                        them once told. Null skips Step F (quantification);
	 *                        peptides then keep phase-1 zero intensity, exactly
	 *                        as every caller of the deprecated overload above
	 *                        already gets.
	 * @param flashLfqThreads thread count passed straight through to FlashLFQ
	 * @param cb              progress sink; may be null
	 */
	public static Result run(DDAMagOpenSearchParameters params, File msfraggerJar, File markerFasta, MagDbItem mag,
			File flashLfqExe, int flashLfqThreads, ProgressCallback cb) throws Exception {
		return run(params, msfraggerJar, markerFasta, mag, flashLfqExe, flashLfqThreads, null, cb);
	}

	/**
	 * Same as the 7-argument {@link #run}, with an additional opt-in
	 * {@code funcDbFile} (MetaPilot's {@code func_def.db}). Null keeps exactly
	 * the 7-argument overload's existing behaviour: functional annotation is
	 * skipped, matching {@code MagFragPipeHeadlessRunner}'s own null-funcDbFile
	 * handling ("FragPipe functional annotation: skipped (func_def.db not
	 * configured)"). This is a brand new overload, not a change to any existing
	 * signature, so it cannot break any existing caller that has not opted in -
	 * same reasoning {@code MagFragPipeHeadlessRunner} used when it added its own
	 * {@code funcDbFile} overload.
	 * <p>
	 * Taxonomy export (Step H below) is <em>not</em> gated on this parameter - it
	 * needs no function-definition database, only the genome-level maps this
	 * workflow already computes - so it now runs unconditionally whenever
	 * post-analysis runs, the same way it always has for FragPipe and pFind.
	 * Before this overload existed, the open-search workflow produced no
	 * {@code taxonomy_analysis/} or {@code functional_annotation/} output at all
	 * (confirmed absent from every call site); every existing caller of the
	 * 7-argument overload above now gets {@code taxonomy_analysis/} as well, and
	 * keeps skipping functional annotation until it is updated to call this
	 * overload with a real {@code funcDbFile}.
	 */
	public static Result run(DDAMagOpenSearchParameters params, File msfraggerJar, File markerFasta, MagDbItem mag,
			File flashLfqExe, int flashLfqThreads, File funcDbFile, ProgressCallback cb) throws Exception {
		return run(params, msfraggerJar, markerFasta, mag, flashLfqExe, flashLfqThreads, funcDbFile, cb, null);
	}

	/**
	 * Same as the 8-argument {@link #run}, with an optional
	 * {@link PartialIdentificationHandler} so the GUI can offer
	 * continue-with-successful when stage-1 evidence is empty for some samples.
	 * Headless callers pass {@code null} and use
	 * {@link #CONTINUE_WITH_SUCCESSFUL_PROPERTY} / CLI
	 * {@code --continue-with-successful-files} instead.
	 */
	public static Result run(DDAMagOpenSearchParameters params, File msfraggerJar, File markerFasta, MagDbItem mag,
			File flashLfqExe, int flashLfqThreads, File funcDbFile, ProgressCallback cb,
			PartialIdentificationHandler partialIdentificationHandler) throws Exception {
		if (params != null && params.isEntrapmentEnabled()) {
			ProgressCallback callback = cb == null ? new ProgressCallback.NoOp() : cb;
			callback.onMessage("Open search: paired-shuffle entrapment evaluation (seed "
					+ params.getEntrapmentSeed() + ", ID digit " + params.getEntrapmentIdDigit() + ")");
			return runEntrapment(params, msfraggerJar, markerFasta, mag, callback,
					new OpenSearchEntrapment.Config(params.getEntrapmentSeed(), params.getEntrapmentIdDigit()));
		}
		return runInternal(params, msfraggerJar, markerFasta, mag, flashLfqExe, flashLfqThreads, funcDbFile, cb, null,
				partialIdentificationHandler);
	}

	/**
	 * Rebuild final tables, taxonomy and functional annotation from an already
	 * finished search. Reads {@code validated_psms.tsv} and
	 * {@code QuantifiedPeptides.tsv}; does not search, localise, classify or
	 * quantify. The genome freeze ({@code MagResultFinalisationService}) and the
	 * current export columns (shift explanation, spectrum count, mzML
	 * identification rate) are applied to those tables.
	 */
	public static Result runPostAnalysisOnly(DDAMagOpenSearchParameters params, MagDbItem mag, File funcDbFile,
			File validatedPsmsOverride, File quantifiedOverride, ProgressCallback cb) throws Exception {
		if (cb == null) {
			cb = new ProgressCallback.NoOp();
		}
		if (params == null) {
			throw new IllegalArgumentException("params is null");
		}
		if (mag == null) {
			throw new IllegalArgumentException("mag is null");
		}
		File userResult = new File(params.getResult());
		if (!userResult.isDirectory()) {
			throw new IOException("Result folder does not exist (pass --out if the parameter file still "
					+ "points at the original HPC path): " + userResult.getAbsolutePath());
		}
		File resultRoot = OpenSearchLayout.engineSearchRoot(userResult);
		File modeRoot = OpenSearchLayout.modeRoot(resultRoot, params.getSearchEngines());
		params.setDbSearchResultFile(modeRoot);
		List<String[]> experiments = params.getExperimentData();
		if (experiments == null || experiments.isEmpty()) {
			throw new IOException("No spectra files configured for the open-search workflow");
		}
		File validatedTsv = validatedPsmsOverride != null ? validatedPsmsOverride
				: OpenSearchLayout.validatedPsms(resultRoot, params.getSearchEngines());
		File quantifiedTsv = quantifiedOverride != null ? quantifiedOverride
				: OpenSearchLayout.quantifiedPeptides(resultRoot, params.getSearchEngines());
		if (quantifiedTsv == null || !quantifiedTsv.isFile() || quantifiedTsv.length() == 0) {
			quantifiedTsv = null;
		}
		cb.onMessage("Resuming post-analysis from " + validatedTsv.getAbsolutePath()
				+ (quantifiedTsv == null ? "" : ", " + quantifiedTsv.getAbsolutePath()));
		OpenSearchValidatedPsmReader.Result loaded = OpenSearchValidatedPsmReader.read(validatedTsv);
		int overlaid = OpenSearchValidatedPsmReader.overlayShifts(loaded.psms,
				OpenSearchLayout.scoredPsms(resultRoot, params.getSearchEngines()));
		if (overlaid > 0) {
			cb.onMessage("Post-analysis: filled shift localisation on " + overlaid + " validated PSMs from "
					+ OpenSearchLayout.scoredPsms(resultRoot, params.getSearchEngines()).getName());
		}
		writeValidatedPsms(loaded.psms, loaded.identificationOnly, validatedTsv);
		MagResultCheckpoints.clearSource(modeRoot);
		exportPostAnalysis(params, mag, resultRoot, modeRoot, experiments, loaded.psms, quantifiedTsv, funcDbFile, cb);
		MagResultCheckpoints.recordSource(modeRoot, validatedTsv, quantifiedTsv);
		cb.onMessage("Open-search post-analysis-only run completed.");
		cb.onComplete(modeRoot.getAbsolutePath());
		return new Result(resultRoot, modeRoot, validatedTsv, OpenSearchLayout.peakVerdicts(resultRoot,
				params.getSearchEngines()), loaded.psms.size(), -1);
	}

	/** Evaluation: add paired shuffles only after genome selection and stop after identification. */
	public static Result runEntrapment(DDAMagOpenSearchParameters params, File msfraggerJar, File markerFasta,
			MagDbItem mag, ProgressCallback cb, OpenSearchEntrapment.Config entrapment) throws Exception {
		if (entrapment == null) throw new IllegalArgumentException("Entrapment configuration is required");
		if (!"Trypsin".equalsIgnoreCase(params.getEnzyme()) && !"Trypsin/P".equalsIgnoreCase(params.getEnzyme())) {
			throw new IllegalArgumentException("Paired entrapment currently supports Trypsin or Trypsin/P");
		}
		params.setClosedStageEnabled(true);
		params.setClosedOpenIntegration(DDAMagOpenSearchParameters.INTEGRATION_JOINT);
		params.setJointCompareArms(true);
		return runInternal(params, msfraggerJar, markerFasta, mag, null, 1, null, cb, entrapment, null);
	}

	private static Result runInternal(DDAMagOpenSearchParameters params, File msfraggerJar, File markerFasta,
			MagDbItem mag, File flashLfqExe, int flashLfqThreads, File funcDbFile, ProgressCallback cb,
			OpenSearchEntrapment.Config entrapment, PartialIdentificationHandler partialIdentificationHandler)
			throws Exception {
		if (cb == null) {
			cb = new ProgressCallback.NoOp();
		}
		if (mag == null) {
			throw new IOException("No MAG database supplied; genome discovery cannot run");
		}
		params.setSearchEngines(DDAMagOpenSearchParameters.ENGINES_MSFRAGGER);
		// Everything lives inside opensearch_MAG_search/ under the user's result
		// folder, as every other MAG workflow does (fragpipe_MAG_search,
		// pfind_MAG_search). Per-sample searches are in search_results/; combined/,
		// final tables, taxonomy, function and report sit beside it. Closed and
		// open MSFragger searches replaced a second engine, so there is no nested
		// msfragger_only/. dbSearchResultFile is this folder because the shared
		// report generator resolves report/ from there.
		File resultRoot = entrapment == null ? OpenSearchLayout.engineSearchRoot(new File(params.getResult()))
				: new File(params.getResult(), OpenSearchEntrapment.RESULT_FOLDER);
		resultRoot.mkdirs();
		OpenSearchLayout.searchResultsRoot(resultRoot).mkdirs();
		File modeRoot = OpenSearchLayout.modeRoot(resultRoot, params.getSearchEngines());
		OpenSearchLayout.combinedDir(resultRoot, params.getSearchEngines()).mkdirs();
		params.setDbSearchResultFile(modeRoot);
		if (!new DDAMagOpenSearchParameterIO().saveAnalysisParameters(modeRoot.getAbsolutePath(), params)) {
			cb.onMessage("WARNING: could not save " + OpenSearchLayout.PARAMETER_FILE + " to " + modeRoot);
		}

		List<String[]> experiments = params.getExperimentData();
		if (experiments == null || experiments.isEmpty()) {
			throw new IOException("No spectra files configured for the open-search workflow");
		}
		// After the parameter file is saved, so it keeps the files as loaded.
		experiments = RawSpectraConversion.prepare(params, cb);

		// --- shared 1:1 target/decoy marker database -----------------------
		// The canonical path from the MAG DB (<current>/Marker4_td.fasta, beside
		// Marker4.fasta), exactly as MagFragPipeHeadlessRunner resolves it. Using
		// the shared location rather than one invented inside the result folder
		// means an existing target-decoy FASTA — and MSFragger's .pepindex built
		// from it — is reused instead of regenerated per run, and both workflows
		// demonstrably search the same file rather than two copies that happen to
		// be built by the same routine.
		File markerTd = mag.getHapFastaTDFile(params.getGenomeSelectorMode());
		if (markerTd == null) {
			// Falls back to the result folder so a MAG DB without a resolved
			// current/ directory still runs, rather than failing outright.
			markerTd = OpenSearchLayout.markerTargetDecoyFasta(resultRoot);
			LOGGER.warn("MAG DB did not supply a target-decoy marker path; falling back to "
					+ markerTd.getAbsolutePath());
		}
		cb.onMessage("Preparing 1:1 target/decoy marker database: " + markerTd.getAbsolutePath());
		if (TargetDecoyFasta.ensure(markerFasta, markerTd)) {
			cb.onMessage("Built target/decoy marker database from " + markerFasta.getName());
		} else {
			cb.onMessage("Reusing existing target/decoy marker database");
		}

		// --- Step A: stage-1 marker search, all samples in one invocation ----
		cb.onProgress(5, 100);
		cb.onMessage("Step A: running " + (params.isRound1ClosedSearch() ? "closed" : "open")
				+ " marker searches (per sample) ...");
		runMarkerSearches(params, msfraggerJar, markerTd, resultRoot, experiments, cb);
		// Precursor refinement before genome selection reads a single delta mass.
		// Stage 1 keeps only the zero family (plus tested UNIMOD shifts in open mode),
		// so a mis-picked precursor is not a small error here - it removes the PSM
		// entirely. See runCrystalC for the measurement. A closed search has no
		// delta-mass distribution to refine: isotope errors are searched
		// directly and every PSM sits within the ppm window.
		if (round1UsesCrystalC(params)) {
			List<File> markerPepXmls = new ArrayList<>(experiments.size());
			for (String[] experiment : experiments) {
				markerPepXmls.add(OpenSearchLayout.markerPepXml(resultRoot, baseName(experiment[0])));
			}
			runCrystalC(params, resultRoot, experiments, markerPepXmls, OpenSearchLayout.MARKER_SEARCH_DIR, markerTd,
					cb);
		} else if (params.isCrystalcEnabled()) {
			cb.onMessage("Crystal-C skipped for the marker search: it is closed");
		}
		// One MSBooster invocation over every sample's pin, before genome
		// selection reads any of them. Batched for the same reason the search is:
		// the model loads once instead of once per sample. Crystal-C does not
		// touch the pin, and the rescorer joins its features back by scan, so the
		// two are independent.
		runStage1MsBooster(params, resultRoot, experiments, cb);

		// --- Step B: zero-peak extraction per sample -----------------------
		cb.onProgress(STEP_B_START_PERCENT, 100);
		cb.onMessage("Step B: selecting genomes and building sample-specific FASTAs ...");
		// The host joins the stage-2 databases only. Entrapment databases pair every real protein
		// with a shuffled MGYG one, which a host protein cannot have, so they stay genomes-only.
		File host = null;
		if (entrapment != null) {
			if (params.isAppendHostDb()) {
				cb.onMessage("Host database not added: entrapment databases hold the selected genomes only");
			}
		} else {
			host = OpenSearchHostDatabase.resolve(params);
			if (host != null) {
				cb.onMessage("Host database appended to each sample-specific database: " + host.getAbsolutePath());
			}
		}
		List<String[]> identificationExperiments = new ArrayList<>();
		List<String> skippedIdentification = new ArrayList<>();
		for (int i = 0; i < experiments.size(); i++) {
			String sample = baseName(experiments.get(i)[0]);
			cb.onMessage("Step B: sample " + (i + 1) + "/" + experiments.size() + " - " + sample);
			if (!extractUnmodified(params, resultRoot, experiments.get(i), mag, host, entrapment != null, cb)) {
				skippedIdentification.add(sample);
			} else {
				identificationExperiments.add(experiments.get(i));
				if (entrapment != null) {
					File targets = refinedTargets(resultRoot, sample);
					OpenSearchEntrapment.prepare(targets, OpenSearchLayout.refinedFasta(resultRoot, sample),
							entrapment, params.sampleDecoyMethod());
					cb.onMessage(sample + ": second-search database includes 1:1 paired shuffled proteins (seed "
							+ entrapment.seed + ", MGYG" + entrapment.idDigit + "...) and "
							+ params.sampleDecoyMethod().label + " decoys of both populations");
				}
			}
			// Advance per sample, as the FragPipe runner does: without this the
			// bar sits still for the whole step and a slow sample is
			// indistinguishable from a hung one.
			double pct = STEP_B_START_PERCENT
					+ (STEP_B_END_PERCENT - STEP_B_START_PERCENT) * (i + 1) / experiments.size();
			cb.onProgress(pct, 100);
		}
		identificationExperiments = applyIdentificationGate(identificationExperiments, skippedIdentification,
				partialIdentificationHandler, cb);

		// --- Step C: stage-2 open search, per sample -----------------------
		// Only samples that built a refined database. Skipped samples still stay
		// in params.metadata for FlashLFQ (Step F), so match-between-runs can
		// quantify them from peptides identified in the rest of the cohort.
		cb.onProgress(45, 100);
		cb.onMessage("Step C: running sample-specific searches (per sample) ...");
		List<File> sampleDbPepXmls = runSampleDbSearches(params, msfraggerJar, resultRoot, identificationExperiments,
				cb);
		// Precursor refinement, before anything reads a delta mass. Everything
		// downstream consumes what this returns; see runCrystalC for why the step
		// exists and what it changes.
		sampleDbPepXmls = runCrystalC(params, resultRoot, identificationExperiments, sampleDbPepXmls,
				OpenSearchLayout.SAMPLE_DB_SEARCH_DIR, null, cb);
		// Sage joins round two only. Round one stays MSFragger-only so both engine
		// configurations select the same genomes and build the same refined
		// databases -- otherwise a difference in the final result could not be
		// attributed to the second engine rather than to a different search space.
		List<File> sageResults = runSageSampleDbSearches(params, resultRoot, identificationExperiments, cb);
		// Closed stage of the closed-then-open identification, on the same refined databases.
		List<File> closedPepXmls = runSampleDbClosedSearches(params, msfraggerJar, resultRoot,
				identificationExperiments, cb);
		cleanupDatabaseIndexes(resultRoot, identificationExperiments, cb);

		// --- Step D: pooled validation -------------------------------------
		cb.onProgress(70, 100);
		String engineMode = params.getSearchEngines();
		File validatedTsv = OpenSearchLayout.validatedPsms(resultRoot, engineMode);
		File scoredTsv = OpenSearchLayout.scoredPsms(resultRoot, engineMode);
		File verdictsTsv = OpenSearchLayout.peakVerdicts(resultRoot, engineMode);
		File peakValidatedTsv = OpenSearchLayout.peakValidatedPsms(resultRoot, engineMode);
		File unexplainedTsv = OpenSearchLayout.unexplainedPsms(resultRoot, engineMode);

		// Staleness again, at the last link in the chain. These five files are
		// derived from the stage-2 pepXMLs, and on 2026-09-11 a run that correctly
		// re-did stage 1 and re-searched every sample then skipped this step
		// because the September 10 outputs still existed - so final_peptides.tsv
		// and razor_proteins.tsv were rebuilt, with fresh timestamps, from PSMs
		// that predated the whole upgrade.
		//
		// Step F (FlashLFQ) and Step G (post-analysis export) did not exist when
		// this gate was written, so it only asked "is validated_psms.tsv itself
		// current". Once they did, that produced a second staleness bug of the
		// same shape: a run where validated_psms.tsv was already current -
		// because nothing upstream of Step D changed - returned here before Step
		// F or G ever ran, silently keeping whatever pre-quantification
		// final_peptides.tsv/razor_proteins.tsv (or none at all) happened to
		// exist. The gate now also requires quantifiedTsv (when a FlashLFQ
		// executable was actually supplied) and the final peptide/protein reports
		// to be present and no older than what feeds them, so "validation is
		// current" only short-circuits the whole method when Steps F and G are
		// ALSO already done - not just Step D.
		List<File> stage2 = new ArrayList<>(sampleDbPepXmls);
		if (sageResults != null) {
			stage2.addAll(sageResults);
		}
		stage2.addAll(closedPepXmls);
		// Switching the closed stage on or off changes no search output the open-only run read,
		// so the switch itself is an input.
		List<File> runSpectra = new ArrayList<>(experiments.size());
		for (String[] experiment : experiments) {
			runSpectra.add(OpenSearchSpectra.readable(new File(experiment[0]), resultRoot));
		}
		stage2.add(writeIdentificationModeStamp(params, resultRoot, engineMode, runSpectra));
		File[] stage2Inputs = stage2.toArray(new File[0]);
		File quantifiedTsv = OpenSearchLayout.quantifiedPeptides(resultRoot, engineMode);
		// Same helpers OpenSearchPostAnalysisTsvExporter.exportFinalTsvs and
		// DDAMagOpenSearchTask (this.final_pep_txt/this.final_pro_txt) both use,
		// rather than duplicating "final_peptides.tsv"/"razor_proteins.tsv" as
		// literals a third time - one shared definition, so this gate can never
		// drift from where those files are actually written or read.
		File finalPeptidesTsv = OpenSearchLayout.finalPeptides(modeRoot);
		File razorProteinsTsv = OpenSearchLayout.razorProteins(modeRoot);
		if (entrapment != null && OpenSearchLayout.isSatisfied(validatedTsv, stage2Inputs)
				&& OpenSearchLayout.isSatisfied(scoredTsv, stage2Inputs)
				&& OpenSearchLayout.isSatisfied(OpenSearchLayout.identificationArms(resultRoot, engineMode), stage2Inputs)
				&& OpenSearchLayout.isSatisfied(OpenSearchLayout.identificationArmPeptides(resultRoot, engineMode), stage2Inputs)) {
			OpenSearchEntrapment.writeSummary(OpenSearchLayout.combinedDir(resultRoot, engineMode), entrapment);
			cb.onMessage("Entrapment identification outputs are current; evaluation summary refreshed");
			cb.onComplete(modeRoot.getAbsolutePath());
			return new Result(resultRoot, modeRoot, validatedTsv, verdictsTsv, countLines(validatedTsv) - 1, -1);
		}
		boolean quantificationExpected = flashLfqExe != null;
		boolean allPresent = OpenSearchLayout.isSatisfied(validatedTsv) && OpenSearchLayout.isSatisfied(scoredTsv)
				&& OpenSearchLayout.isSatisfied(verdictsTsv) && OpenSearchLayout.isSatisfied(peakValidatedTsv)
				&& OpenSearchLayout.isSatisfied(unexplainedTsv)
				&& (!quantificationExpected || OpenSearchLayout.isSatisfied(quantifiedTsv))
				&& OpenSearchLayout.isSatisfied(finalPeptidesTsv) && OpenSearchLayout.isSatisfied(razorProteinsTsv);
		if (allPresent) {
			boolean current = OpenSearchLayout.isSatisfied(validatedTsv, stage2Inputs)
					&& OpenSearchLayout.isSatisfied(scoredTsv, stage2Inputs)
					&& OpenSearchLayout.isSatisfied(verdictsTsv, stage2Inputs)
					&& OpenSearchLayout.isSatisfied(peakValidatedTsv, stage2Inputs)
					&& OpenSearchLayout.isSatisfied(unexplainedTsv, stage2Inputs)
					&& (!quantificationExpected || OpenSearchLayout.isSatisfied(quantifiedTsv, validatedTsv))
					&& OpenSearchLayout.isSatisfied(finalPeptidesTsv, quantifiedTsv, validatedTsv)
					&& OpenSearchLayout.isSatisfied(razorProteinsTsv, quantifiedTsv, validatedTsv)
					&& MagResultCheckpoints.sameSource(modeRoot, validatedTsv,
							quantificationExpected ? quantifiedTsv : null);
			// The pangenome variant check needs this run's SubstitutionReassigner
			// result, which exists only in memory during Step D. When the add-on
			// would run and its outputs are missing or older than what they were
			// derived from, skipping here would leave it never written.
			boolean pangenomeCurrent = !current
					|| pangenomeVariantOutputsCurrent(params, mag, resultRoot, modeRoot, validatedTsv);
			if (current && pangenomeCurrent) {
				cb.onMessage("Validated PSMs, quantification and final reports already exist and are current, "
						+ "skipping straight to completion: " + validatedTsv.getName());
				quantifySubstitutionPsms(params, flashLfqExe, flashLfqThreads, resultRoot, engineMode,
						readableSpectraDir(experiments, resultRoot), cb);
				return new Result(resultRoot, modeRoot, validatedTsv, verdictsTsv, countLines(validatedTsv) - 1,
						-1);
			}
			if (current) {
				cb.onMessage("Pangenome variant check outputs are missing or older than validated_psms.tsv; "
						+ "re-running validation so they can be written (quantification and export are kept "
						+ "when validated_psms.tsv comes out unchanged)");
			} else {
				File newest = OpenSearchLayout.newestInput(stage2Inputs);
				cb.onMessage("Validated PSMs, quantification or final reports are missing, older than "
						+ (newest == null ? "the stage-2 searches" : newest.getName())
						+ " or not recorded as produced from " + validatedTsv.getName()
						+ "; re-running validation, quantification and export rather than reporting a stale result");
			}
		}

		cb.onMessage("Step D: delta-mass fitting, classification and per-modification validation ...");
		// Step D is deterministic, so a rerun forced only by missing report-only
		// outputs rewrites validated_psms.tsv with the same bytes. Keeping its
		// timestamp then lets Steps F-I see unchanged input and keep their outputs.
		MagResultCheckpoints.UnchangedRewrite validatedBefore = MagResultCheckpoints.UnchangedRewrite
				.before(validatedTsv, stage2Inputs);
		AtomicReference<SubstitutionReassigner.Result> substitutions = new AtomicReference<>();
		List<File> spectraFiles = new ArrayList<>(experiments.size());
		for (String[] experiment : experiments) {
			spectraFiles.add(OpenSearchSpectra.readable(new File(experiment[0]), resultRoot));
		}
		ModificationValidator.Result validation = validatePepXmls(params, sampleDbPepXmls, sageResults, closedPepXmls,
				spectraFiles, resultRoot, cb, substitutions);
		if (entrapment != null) {
			OpenSearchEntrapment.writeSummary(OpenSearchLayout.combinedDir(resultRoot, engineMode), entrapment);
			cb.onMessage("Entrapment evaluation completed at identification; controls retained in combined tables");
			cb.onProgress(100, 100);
			cb.onComplete(modeRoot.getAbsolutePath());
			return new Result(resultRoot, modeRoot, validatedTsv, verdictsTsv, validation.validated.size(),
					validation.acceptedPeakCount());
		}
		if (validatedBefore.restoreIfUnchanged()) {
			cb.onMessage("Step D: " + validatedTsv.getName() + " is unchanged; kept its timestamp so the "
					+ "quantification and reports derived from it stay current");
		} else if (validatedBefore.remembered()) {
			cb.onMessage("Step D: " + validatedTsv.getName() + " changed; quantification and reports derived "
					+ "from it will be regenerated");
		}

		cb.onProgress(88, 100);
		cb.onMessage("Step F: FlashLFQ quantification ...");
		quantifiedTsv = runFlashLfqQuantification(params, flashLfqExe, flashLfqThreads, resultRoot,
				OpenSearchLayout.validatedPsms(resultRoot, engineMode), readableSpectraDir(experiments, resultRoot), engineMode,
				cb);

		boolean functionExpected = funcDbFile != null && funcDbFile.isFile();
		File stalePostAnalysis = MagResultCheckpoints.stalePostAnalysisOutput(modeRoot, functionExpected,
				validatedTsv, quantifiedTsv);
		boolean sameSource = MagResultCheckpoints.sameSource(modeRoot, validatedTsv, quantifiedTsv);
		File qpxPsms = new File(modeRoot, QpxPsmParquetWriter.DEFAULT_FILE_NAME);
		if (stalePostAnalysis == null && sameSource) {
			cb.onProgress(90, 100);
			cb.onMessage("Steps G-I: final tables, taxonomy" + (functionExpected ? " and functional annotation" : "")
					+ " were produced from " + validatedTsv.getName() + " and are newer than it"
					+ (quantifiedTsv == null ? "" : " and " + quantifiedTsv.getName()) + ", skipping");
			if (!OpenSearchLayout.isSatisfied(qpxPsms, validatedTsv)) {
				List<String> sampleNames = new ArrayList<>();
				for (String[] experiment : experiments) {
					sampleNames.add(baseName(experiment[0]));
				}
				OpenSearchPostAnalysisTsvExporter.writeQpxPsms(validation.validated, sampleNames, resultRoot,
						modeRoot, cb);
			}
		} else {
			if (stalePostAnalysis != null) {
				cb.onMessage("Steps G-I: " + stalePostAnalysis.getName() + " is missing or older than its inputs; "
						+ "regenerating final tables, taxonomy and functional annotation");
			} else {
				cb.onMessage("Steps G-I: the final tables, taxonomy and functional annotation in "
						+ modeRoot.getName() + " were not recorded as produced from " + validatedTsv.getName()
						+ " (an interrupted export, or they predate this record); regenerating");
			}
			// A parquet left from an earlier run would otherwise outlive a failed export.
			MagResultCheckpoints.clearSource(modeRoot);
			Files.deleteIfExists(qpxPsms.toPath());
			exportPostAnalysis(params, mag, resultRoot, modeRoot, experiments, validation.validated, quantifiedTsv,
					funcDbFile, cb);
			MagResultCheckpoints.recordSource(modeRoot, validatedTsv, quantifiedTsv);
		}

		runPangenomeVariantCheck(params, mag, resultRoot, modeRoot, validation, substitutions.get(), cb);
		quantifySubstitutionPsms(params, flashLfqExe, flashLfqThreads, resultRoot, engineMode,
				readableSpectraDir(experiments, resultRoot), cb);

		LOGGER.info("Open-search workflow complete: " + validation.validated.size() + " validated PSMs across "
				+ validation.acceptedPeakCount() + " accepted delta-mass groups");
		cb.onProgress(100, 100);
		cb.onMessage("MAG+MSFragger open-search headless run completed.");
		cb.onComplete(modeRoot.getAbsolutePath());

		return new Result(resultRoot, modeRoot, validatedTsv, verdictsTsv, validation.validated.size(),
				validation.acceptedPeakCount());
	}

	/**
	 * Steps G-I: post-analysis export, taxonomy and functional annotation. One
	 * unit because taxonomy and function read the in-memory result of the
	 * export.
	 *
	 * @param resultRoot the engine search root ({@code opensearch_MAG_search/})
	 * @param modeRoot   same folder; kept as a separate argument for call sites
	 *                   that still distinguish searches from deliverables
	 */
	private static void exportPostAnalysis(DDAMagOpenSearchParameters params, MagDbItem mag, File resultRoot,
			File modeRoot, List<String[]> experiments, List<OpenSearchPsm> validated, File quantifiedTsv,
			File funcDbFile, ProgressCallback cb) throws Exception {
		cb.onProgress(90, 100);
		cb.onMessage("Step G: post-analysis + TSV export ...");
		// Hands the validated PSMs to the same generic services every DDA MAG
		// engine uses (genome scoring, protein inference, razor assignment) and
		// writes the standard deliverables.
		List<String> sampleNames = new ArrayList<>();
		for (String[] experiment : experiments) {
			sampleNames.add(baseName(experiment[0]));
		}
		// Identification v2: RESET's estimating decoys count twice (C2), one-peptide protein filter (C5).
		GenomeEvidenceTest.Result genomeEvidence = OpenSearchGenomeEvidence.run(resultRoot, params.getSearchEngines(),
				sampleNames, cb, params.isIdentificationV2() ? 2.0 : 1.0, params.isIdentificationV2(),
				mag.getIdentifier());
		OpenSearchPostAnalysisTsvExporter.Result post = OpenSearchPostAnalysisTsvExporter.exportFinalTsvs(mag,
				validated, sampleNames, modeRoot, resultRoot, spectraSummaryDir(experiments, resultRoot, cb),
				"DDA MAG open search", quantifiedTsv, genomeEvidence, cb);
		cb.onMessage("Post-analysis: " + post.peptideCount + " peptides, " + post.proteinCount + " razor proteins, "
				+ post.genomeScoreMap.size() + " genomes");
		if (post.peptidesWithoutConfidentRazorProtein > 0) {
			cb.onMessage("Post-analysis: " + post.peptidesWithoutConfidentRazorProtein + " of " + post.peptideCount
					+ " peptides have no confident genome/razor-protein attribution (every candidate genome is at "
					+ "q_value > " + DdaMagProteinInferenceService.DEFAULT_GENOME_Q_VALUE_THRESHOLD + ") - still "
					+ "recorded in excluded_evidence/peptides.tsv and excluded from every "
					+ "genome- and protein-level statistic");
		}

		cb.onMessage("Step H: taxonomy report export ...");
		exportTaxonomyTsvs(mag, modeRoot, post, cb);

		if (funcDbFile != null && funcDbFile.isFile()) {
			cb.onMessage("Step I: functional annotation report export ...");
			MagFunctionalAnnotationService.Request funcRequest = new MagFunctionalAnnotationService.Request();
			funcRequest.resultDir = modeRoot;
			funcRequest.magDb = mag;
			funcRequest.microDbCatalogueId = params.getMagDbCatalogueId();
			funcRequest.functionDefinitionDb = funcDbFile;
			funcRequest.analysisType = MetaConstants.openSearchMAG;
			funcRequest.quantitativeMode = params.getQuantitativeMode();
			funcRequest.intensityTitles = post.intensityTitles;
			funcRequest.metaProteins = post.metaProteins;
			MagFunctionalAnnotationService.run(funcRequest);
		} else {
			cb.onMessage("Open search functional annotation: skipped (func_def.db not configured)");
		}
	}

	/**
	 * Pangenome variant check (UHGP-100), after taxonomy so the DETECTED species
	 * scope can use {@code genome_final_selected.tsv}. Report-only: writes the
	 * {@code variant_*.tsv} files in {@code sequence_variants/} and
	 * changes nothing else. Off, or with
	 * no resource present, it writes nothing; any failure is logged and the run
	 * still completes.
	 */
	private static void runPangenomeVariantCheck(DDAMagOpenSearchParameters params, MagDbItem mag, File resultRoot,
			File modeRoot, ModificationValidator.Result validation, SubstitutionReassigner.Result substitutions,
			ProgressCallback cb) {
		try {
			if (params.getPangenomeVariantSettings()
					.getMode() == uk.ac.quadram.metapilot.core.variant.PangenomeVariantSettings.Mode.OFF) {
				return;
			}
			String engineMode = params.getSearchEngines();
			OpenSearchVariantCandidateSource.Result source = OpenSearchVariantCandidateSource.collect(substitutions,
					validation.diagnostic, params.getSubstitutionMassToleranceDa());
			PangenomeVariantCheck.Request request = new PangenomeVariantCheck.Request();
			request.engine = VariantCandidate.Engine.OPEN_SEARCH;
			request.candidates = source.candidates;
			request.identifiedBasePeptides = source.basePeptides;
			request.searchedFastas = refinedFastas(params, resultRoot);
			request.decoyPrefix = TargetDecoyFasta.DECOY_TAG;
			request.settings = params.getPangenomeVariantSettings();
			request.magVersionDir = mag == null ? null : mag.getCurrentFile();
			request.taxonomyDir = new File(modeRoot, "taxonomy_analysis");
			request.eggnogDir = mag == null ? null : mag.getOriginalFuncFolder();
			request.outputDir = OpenSearchLayout.sequenceVariantsDir(resultRoot, engineMode);
			request.undecidableToleranceDa = params.getSubstitutionMassToleranceDa();
			request.fixedCysteineModDa = SubstitutionMassTable.FIXED_CYSTEINE_MOD;
			VariantFormRanking.Inputs ranking = new VariantFormRanking.Inputs();
			ranking.apply(params.getVariantFormRankingSettings());
			ranking.diann = params.getDiannExecutable();
			ranking.threads = Math.max(1, params.getThreadCount());
			ranking.substitutionToleranceDa = params.getSubstitutionMassToleranceDa();
			List<String[]> experiments = params.getExperimentData();
			if (experiments != null) {
				for (String[] experiment : experiments) {
					ranking.spectraFiles.add(OpenSearchSpectra.readable(new File(experiment[0]), resultRoot));
				}
			}
			request.ranking = ranking;
			request.info = msg -> {
				LOGGER.info(msg);
				cb.onMessage(msg);
			};
			request.warn = msg -> {
				LOGGER.warn(msg);
				cb.onMessage("WARNING: " + msg);
			};
			if (PangenomeVariantCheck.wouldRun(request.settings, request.magVersionDir)) {
				cb.onMessage("Pangenome variant check: " + source.describe());
			}
			PangenomeVariantCheck.run(request);
		} catch (RuntimeException e) {
			LOGGER.warn("Pangenome variant check failed; workflow results are unaffected", e);
			cb.onMessage("WARNING: Pangenome variant check failed; workflow results are unaffected: " + e);
		}
	}

	/**
	 * False only when the pangenome variant check would run and its outputs are
	 * missing or older than the validated PSMs or the genome selection they
	 * were derived from.
	 */
	private static boolean pangenomeVariantOutputsCurrent(DDAMagOpenSearchParameters params, MagDbItem mag,
			File resultRoot, File modeRoot, File validatedTsv) {
		if (!PangenomeVariantCheck.wouldRun(params.getPangenomeVariantSettings(),
				mag == null ? null : mag.getCurrentFile())) {
			return true;
		}
		String engineMode = params.getSearchEngines();
		OpenSearchLayout.variantColumns(resultRoot, engineMode);
		File candidates = OpenSearchLayout.variantCandidates(resultRoot, engineMode);
		File sites = OpenSearchLayout.variantSites(resultRoot, engineMode);
		File summary = OpenSearchLayout.variantPangenomeSummary(resultRoot, engineMode);
		File species = OpenSearchLayout.variantSpeciesSummary(resultRoot, engineMode);
		File selected = new File(new File(modeRoot, "taxonomy_analysis"), "genome_final_selected.tsv");
		return OpenSearchLayout.isSatisfied(candidates, validatedTsv)
				&& OpenSearchLayout.isSatisfied(sites, candidates, validatedTsv)
				&& OpenSearchLayout.isSatisfied(summary, candidates, validatedTsv, selected)
				&& OpenSearchLayout.isSatisfied(species, candidates, validatedTsv, selected);
	}

	/**
	 * Step H: shared taxonomy analysis, mirroring
	 * {@code MagFragPipeHeadlessRunner#exportTaxonomyTsvs} field-for-field - same
	 * source ({@link OpenSearchPostAnalysisTsvExporter.Result}, itself built from
	 * the same {@code DdaMagPostAnalysisService.Result} FragPipe's exporter
	 * reads), same target ({@link MagTaxonomyAnalysisService}), so
	 * {@code taxonomy_analysis/} looks the same regardless of which engine
	 * produced it. Unlike functional annotation, this needs no external
	 * function-definition database and so is not gated on one - only on
	 * post-analysis having run at all, which by the time this is called it has.
	 */
	private static void exportTaxonomyTsvs(MagDbItem mag, File modeRoot, OpenSearchPostAnalysisTsvExporter.Result post,
			ProgressCallback cb) throws IOException {
		File taxDir = new File(modeRoot, "taxonomy_analysis");
		taxDir.mkdirs();

		File gtdbNcbi = new File("Resources", "gtdb_ncbi");
		if (!gtdbNcbi.exists()) {
			gtdbNcbi = null;
		}
		MagTaxonomyAnalysisService.Request request = new MagTaxonomyAnalysisService.Request();
		request.taxonomyAnalysisDir = taxDir;
		request.magDb = mag;
		request.gtdbNcbiMappingTsv = gtdbNcbi;
		request.generateFinalGenomeRanking = true;
		request.genomeIntensityMap = post.genomeIntensityMap;
		request.genomeUniqueIntensityMap = post.genomeUniqueIntensityMap;
		request.genomeScoreMap = post.genomeScoreMap;
		request.genomeStats = post.genomeDiaStats;
		request.genomeEvidence = post.genomeEvidence;
		request.psmCountMap = post.psmCountMap;
		request.totalPepCountMap = post.totalPepCountMap;
		request.razorPepCountMap = post.razorPepCountMap;
		request.uniquePepCountMap = post.uniquePepCountMap;
		request.totalPepWithILMap = post.totalPepWithILMap;
		request.razorPepWithILMap = post.razorPepWithILMap;
		request.uniquePepWithILMap = post.uniquePepWithILMap;
		request.intensityTitles = post.intensityTitles;
		request.progressCallback = cb;
		MagTaxonomyAnalysisService.run(request);
	}

	/**
	 * Step F: quantifies {@code validated_psms.tsv} with FlashLFQ, a single
	 * combined run over every sample together — mirroring
	 * {@code DDAMagPfindTask.lfQuant()}'s
	 * {@code flashLfqTask.run(magPfindParams.getMetadata(), msms.getAbsolutePath(),
	 * rawDir.getAbsolutePath(), true, true)} exactly, down to {@code mbr}/{@code nor}
	 * both {@code true}. Anything else — quantifying per sample, or skipping
	 * FlashLFQ's own combined-run cross-sample reconciliation — would be a new,
	 * unproven mechanism rather than the one pFind's workflow already
	 * establishes.
	 * <p>
	 * Only a missing {@code flashLfqExe} degrades gracefully (returns null, so
	 * {@link OpenSearchPostAnalysisTsvExporter#exportFinalTsvs} keeps phase-1's
	 * zero-intensity peptides) — that is what a caller not supplying one is
	 * asking for, the {@linkplain #run(DDAMagOpenSearchParameters, File, File,
	 * MagDbItem, ProgressCallback) deprecated identify-only overload}'s whole
	 * point. Once FlashLFQ is actually invoked, a failure throws rather than
	 * silently degrading, matching {@code DDAMagPfindTask.lfQuant()}: it checks
	 * {@code lfqQuan} and returns {@code false} (failing the whole task) rather
	 * than proceeding with unquantified data. A "full test" run should surface a
	 * quantification failure loudly, not produce a report that looks complete
	 * but silently has every peptide at zero intensity.
	 */
	private static File runFlashLfqQuantification(DDAMagOpenSearchParameters params, File flashLfqExe,
			int flashLfqThreads, File resultRoot, File validatedTsv, File spectraDir, String engineMode,
			ProgressCallback cb) throws IOException {
		File quantifiedTsv = OpenSearchLayout.quantifiedPeptides(resultRoot, engineMode);
		if (flashLfqExe == null || !flashLfqExe.isFile()) {
			cb.onMessage("Step F: no FlashLFQ executable supplied; skipping quantification, "
					+ "peptides will keep phase-1 zero intensity");
			return null;
		}
		if (!OpenSearchLayout.isSatisfied(validatedTsv)) {
			throw new IOException(
					"Step F: " + validatedTsv.getAbsolutePath() + " is missing or empty; cannot quantify");
		}
		if (spectraDir == null || !spectraDir.isDirectory()) {
			throw new IOException("Step F: could not resolve the raw/mzML directory (" + spectraDir
					+ "); cannot quantify");
		}
		if (OpenSearchLayout.isSatisfied(quantifiedTsv, validatedTsv)) {
			cb.onMessage("Step F: " + quantifiedTsv.getName() + " already exists and is current, skipping FlashLFQ");
			return quantifiedTsv;
		}
		MetaData metadata = params.getMetadata();
		if (metadata == null || metadata.getRawFileCount() == 0) {
			throw new IOException("Step F: workflow parameters carry no raw-file metadata; cannot quantify");
		}
		cb.onMessage("Step F: running FlashLFQ over " + validatedTsv.getName() + " and " + spectraDir);
		FlashLfqTask task = new FlashLfqTask(flashLfqExe.getAbsolutePath(), flashLfqThreads);
		boolean ok = task.run(metadata, validatedTsv.getAbsolutePath(), spectraDir.getAbsolutePath(), true, true);
		if (!ok) {
			throw new IOException("Step F: FlashLFQ quantification failed for " + validatedTsv.getAbsolutePath()
					+ " - see logging.log for FlashLFQ's own output");
		}
		if (!OpenSearchLayout.isSatisfied(quantifiedTsv)) {
			throw new IOException("Step F: FlashLFQ reported success but " + quantifiedTsv.getAbsolutePath()
					+ " is missing or empty");
		}
		cb.onMessage("Step F: FlashLFQ quantification complete: " + quantifiedTsv);
		return quantifiedTsv;
	}

	/**
	 * Runs only the pooled validation stage against existing MSFragger pepXML
	 * files. This is the reproducible entry point for evaluation/upgrade loops:
	 * it never launches a search and never modifies the source pepXML files.
	 * Outputs are written below {@code resultRoot/combined} using the normal
	 * workflow filenames.
	 */
	public static ModificationValidator.Result validatePepXmls(DDAMagOpenSearchParameters params,
			List<File> pepXmlFiles, File resultRoot, ProgressCallback cb) throws Exception {
		return validatePepXmls(params, pepXmlFiles, Collections.emptyList(), resultRoot, cb);
	}

	/**
	 * @param sageResults per-sample {@code results.sage.tsv} files, or empty for a
	 *                    MSFragger-only run. Sage candidates are appended to the
	 *                    same pool and go through every downstream step
	 *                    identically; only the classifier is engine-aware, because
	 *                    Sage cannot supply the delta-mass localization features.
	 */
	public static ModificationValidator.Result validatePepXmls(DDAMagOpenSearchParameters params,
			List<File> pepXmlFiles, List<File> sageResults, File resultRoot, ProgressCallback cb) throws Exception {
		return validatePepXmls(params, pepXmlFiles, sageResults, resultRoot, cb, null);
	}

	/**
	 * @param substitutionsOut when non-null, receives this validation's
	 *                         {@link SubstitutionReassigner} result, for the
	 *                         report-only pangenome variant check
	 */
	public static ModificationValidator.Result validatePepXmls(DDAMagOpenSearchParameters params,
			List<File> pepXmlFiles, List<File> sageResults, File resultRoot, ProgressCallback cb,
			AtomicReference<SubstitutionReassigner.Result> substitutionsOut) throws Exception {
		return validatePepXmls(params, pepXmlFiles, sageResults, null, resultRoot, cb, substitutionsOut);
	}

	/**
	 * @param closedPepXmls the closed searches of the closed-then-open identification; it runs when
	 *                      these are supplied and the closed stage is enabled, otherwise the
	 *                      open-only workflow does
	 */
	public static ModificationValidator.Result validatePepXmls(DDAMagOpenSearchParameters params,
			List<File> pepXmlFiles, List<File> sageResults, List<File> closedPepXmls, File resultRoot,
			ProgressCallback cb, AtomicReference<SubstitutionReassigner.Result> substitutionsOut) throws Exception {
		return validatePepXmls(params, pepXmlFiles, sageResults, closedPepXmls, null, resultRoot, cb, substitutionsOut);
	}

	/**
	 * @param spectraFiles the runs' spectrum files; the joint identification reads the MS2 peaks
	 *                     of {@code <sample>.mzML} for its predicted-spectrum features. Null or
	 *                     empty turns those features off.
	 */
	public static ModificationValidator.Result validatePepXmls(DDAMagOpenSearchParameters params,
			List<File> pepXmlFiles, List<File> sageResults, List<File> closedPepXmls, List<File> spectraFiles,
			File resultRoot, ProgressCallback cb, AtomicReference<SubstitutionReassigner.Result> substitutionsOut)
			throws Exception {
		if (params == null) {
			params = new DDAMagOpenSearchParameters();
		}
		if (pepXmlFiles == null || pepXmlFiles.isEmpty()) {
			throw new IOException("No pepXML files supplied for pooled validation");
		}
		if (resultRoot == null) {
			throw new IOException("No validation output folder supplied");
		}
		if (cb == null) {
			cb = new ProgressCallback.NoOp();
		}
		for (File pepXml : pepXmlFiles) {
			if (pepXml == null || !pepXml.isFile() || pepXml.length() == 0) {
				throw new IOException("Missing or empty pepXML: " + pepXml);
			}
		}
		resultRoot.mkdirs();
		OpenSearchLayout.combinedDir(resultRoot, params.getSearchEngines()).mkdirs();

		OpenSearchPepxmlReader reader = new OpenSearchPepxmlReader(TargetDecoyFasta.DECOY_TAG);
		final ProgressCallback validationCb = cb;
		List<OpenSearchPsm> psms = reader.read(pepXmlFiles, validationCb::onMessage);
		if (psms.isEmpty()) {
			throw new IOException("Supplied pepXML files produced no PSMs; cannot validate");
		}
		int msfraggerCount = psms.size();
		if (sageResults != null && !sageResults.isEmpty()) {
			SageOpenSearchReader sageReader = new SageOpenSearchReader(TargetDecoyFasta.DECOY_TAG);
			for (File sageTsv : sageResults) {
				String sample = sageTsv.getParentFile() == null ? sageTsv.getName()
						: sageTsv.getParentFile().getParentFile().getName();
				psms.addAll(sageReader.read(sageTsv, sample, validationCb::onMessage));
			}
			sageReader.reportUnrecognisedColumns(validationCb::onMessage);
			cb.onMessage("Pooled candidates: " + msfraggerCount + " MSFragger + "
					+ (psms.size() - msfraggerCount) + " Sage = " + psms.size());
		}
		// Missing localization evidence degrades the run; it does not invalidate
		// it. The classifier keeps its search, delta-mass and peak-context
		// features and still produces a scored, FDR-controlled result — just a
		// weaker one. Aborting here yields nothing to inspect and no way to tell
		// how much the localization features were actually contributing, so warn
		// unmissably and continue. The score-name list logged just above is what
		// distinguishes "MSFragger did not emit it" from "emitted under a name
		// this reader does not map".
		if (!OpenSearchPepxmlReader.hasLocalizationEvidence(psms)) {
			String warning = "No delta-mass localization evidence in the pepXML. The classifier will run WITHOUT "
					+ "its localization features (score/ions with and without the delta mass, best localized "
					+ "position). Results are provisional: PSM scores and FDR remain valid, but modification "
					+ "site-level confidence is not being assessed. Check localize_delta_mass = 1 in "
					+ OpenSearchLayout.SEARCH_PARAMS_FILE + " beside the pepXML, and inspect the ptm_result "
					+ "attribute names listed above.";
			cb.onMessage("WARNING: " + warning);
			LOGGER.warn(warning);
		}

		List<OpenSearchPsm> closedPsms = null;
		if (params.isClosedStageEnabled() && closedPepXmls != null && !closedPepXmls.isEmpty()) {
			for (File pepXml : closedPepXmls) {
				if (pepXml == null || !pepXml.isFile() || pepXml.length() == 0) {
					throw new IOException("Missing or empty closed-search pepXML: " + pepXml);
				}
			}
			OpenSearchPepxmlReader closedReader = new OpenSearchPepxmlReader(TargetDecoyFasta.DECOY_TAG);
			closedReader.setClosedSearch(true);
			// Every closed candidate the search reported (top closed_stage_topn) is read, in every mode;
			// the cascade's closed stage still keeps each spectrum's best-ranked one.
			closedReader.setMaxHitRank(params.getClosedStageTopN());
			closedPsms = closedReader.read(closedPepXmls, validationCb::onMessage);
		}
		ModificationValidator.Result validation = validate(params, psms, closedPsms, spectraFiles, resultRoot, cb,
				substitutionsOut);
		// Must run before anything below reads getCanonicalModifiedSequence() on
		// any of these PSMs (writeValidatedPsms's "Full Sequence" column, then
		// later aggregatePeptides and OpenSearchMetaPsm in exportFinalTsvs): it
		// decides, once per peptide+delta-mass group rather than per spectrum,
		// where - if anywhere - the [+mass] bracket sits, so localization
		// ambiguity across spectra can never split one peptide into multiple
		// validated_psms.tsv/QuantifiedPeptides.tsv rows. See
		// OpenSearchPsm.getCanonicalModifiedSequence()'s javadoc.
		OpenSearchPostAnalysisTsvExporter.resolveCanonicalModifiedSequences(reportedPsms(validation));
		// Which walk retained a PSM is a fact about the validation result, not about
		// its peak id: a PSM demoted from a rejected peak keeps peakId >= 0 and is
		// still an identification-only retention.
		Set<OpenSearchPsm> identificationOnly = Collections.newSetFromMap(new IdentityHashMap<>());
		identificationOnly.addAll(validation.unexplainedIdentifications);
		// Each sample searches its own round-2 refined database, so a single PSM's
		// own alternative_protein list is only what THAT sample's database happened
		// to contain - measured at 9.90% row-to-row disagreement for a shared
		// backbone between just two samples. This is NOT reconciled here, matching
		// the pFind workflow exactly: DDAMagPfindTask.combineResults() writes each
		// sample's PFindPSM proteins straight through, per row, with no cross-sample
		// union either. Consistency is FlashLFQ's job, not this file's: its own
		// output is a peptide-level table with exactly one "Protein Groups" value
		// per peptide (FlashLfqQuanPepReader.parse()), aggregated by FlashLFQ across
		// every PSM for that peptide across every sample - and that output, not this
		// file's per-row Protein Accession, is what DdaMagProteinInferenceService
		// actually reads for protein inference. Reconciling it here too would be
		// redundant with what FlashLFQ already guarantees downstream.
		writeValidatedPsms(validation.validated, identificationOnly,
				OpenSearchLayout.validatedPsms(resultRoot, params.getSearchEngines()));
		// Reporting only: groups the substitution PSMs whose candidate is absent
		// from every searched database into protein sites and classes them by
		// recurrence, codon reachability and whether the reference peptide is
		// also seen. Changes no PSM, q-value or quantification input.
		SequenceVariantCaller.Result variants = callSequenceVariants(params, resultRoot, validation.diagnostic,
				unmodifiedIdentifications(validation), cb);
		writeSubstitutionForms(validation.diagnostic, variants,
				OpenSearchLayout.substitutionFormsOutput(resultRoot, params.getSearchEngines()));
		Map<OpenSearchPsm, String> variantIds = new IdentityHashMap<>();
		if (variants != null) {
			variants.byPsm.forEach((psm, v) -> variantIds.put(psm, v.id));
		}
		SubstitutionQuantification.writePsms(validation.diagnostic, validation.validated, variantIds,
				OpenSearchLayout.substitutionQuantPsms(resultRoot, params.getSearchEngines()));
		writeValidatedPsms(validation.peakValidated, identificationOnly,
				OpenSearchLayout.peakValidatedPsms(resultRoot, params.getSearchEngines()));
		writeValidatedPsms(validation.unexplainedIdentifications, identificationOnly,
				OpenSearchLayout.unexplainedPsms(resultRoot, params.getSearchEngines()));
		// The same PSMs per peptide and per genome: the levels the main result is read at.
		int setAsidePeptides = SetAsideReport.write(validation.unexplainedIdentifications, validation.validated,
				OpenSearchLayout.combinedDir(resultRoot, params.getSearchEngines()), TargetDecoyFasta.DECOY_TAG);
		cb.onMessage("Set aside: " + validation.unexplainedIdentifications.size() + " PSMs that are not identifications, "
				+ setAsidePeptides + " peptide rows in " + SetAsideReport.PEPTIDES_FILE + "; per genome in "
				+ SetAsideReport.GENOMES_FILE);
		writeVerdicts(validation, OpenSearchLayout.peakVerdicts(resultRoot, params.getSearchEngines()));
		return validation;
	}

	/** Validated plus diagnostic PSMs, each once: every PSM a written table names. */
	private static List<OpenSearchPsm> reportedPsms(ModificationValidator.Result validation) {
		if (validation.diagnostic == validation.validated) {
			return validation.validated;
		}
		Set<OpenSearchPsm> seen = Collections.newSetFromMap(new IdentityHashMap<>());
		List<OpenSearchPsm> all = new ArrayList<>();
		for (List<OpenSearchPsm> list : Arrays.asList(validation.validated, validation.diagnostic)) {
			for (OpenSearchPsm psm : list) {
				if (seen.add(psm)) {
					all.add(psm);
				}
			}
		}
		return all;
	}

	private static void runMarkerSearches(DDAMagOpenSearchParameters params, File msfraggerJar, File markerTd,
			File resultRoot, List<String[]> experiments, ProgressCallback cb) throws Exception {
		MSFraggerHeadlessRunner.OpenSearchSettings settings = round1Settings(params, markerTd);
		// .raw/.d inputs: the marker search also writes the mzML every later step reads
		// (OpenSearchSpectra). A run whose mzML is missing searches again to write it.
		List<File> inputs = new ArrayList<>(experiments.size());
		for (String[] experiment : experiments) {
			inputs.add(new File(experiment[0]));
		}
		boolean convert = OpenSearchSpectra.needsConversion(inputs);
		settings.writeUncalibratedMzml = convert;
		// Compared with the parameters each existing search was run with, so a
		// change of round-1 settings (closed <-> open, tolerance, top N) re-runs
		// the search instead of reusing a pepXML from the other configuration.
		String wantedParams = MSFraggerHeadlessRunner.previewOpenSearchParamsText(settings);
		List<File> pending = new ArrayList<>();
		List<String> pendingSamples = new ArrayList<>();
		for (String[] experiment : experiments) {
			String sample = baseName(experiment[0]);
			File sentinel = OpenSearchLayout.markerPepXml(resultRoot, sample);
			File existingParams = OpenSearchLayout.markerParams(resultRoot, sample);
			boolean paramsUnchanged = existingParams.isFile()
					&& wantedParams.trim().equals(readFileText(existingParams).trim());
			boolean spectraReadable = OpenSearchSpectra.isReadable(new File(experiment[0]), resultRoot);
			if (paramsUnchanged && OpenSearchLayout.isSatisfied(sentinel) && spectraReadable) {
				cb.onMessage("Marker search already exists under the same configuration, skipping: " + sample);
				continue;
			}
			if (OpenSearchLayout.isSatisfied(sentinel)) {
				cb.onMessage(paramsUnchanged
						? "Marker search re-running to write the mzML of a vendor-format input: " + sample
						: "Marker search configuration changed, re-running: " + sample);
			}
			OpenSearchLayout.markerDir(resultRoot, sample).mkdirs();
			pending.add(new File(experiment[0]));
			pendingSamples.add(sample);
		}
		if (pending.isEmpty()) {
			return;
		}

		// Parameters are written into the FIRST pending sample's marker folder and
		// copied into the others, so every search folder holds both its pepXML and
		// the parameters that produced it and can be checked without looking
		// elsewhere in the tree. Round one is a single batch invocation, so one
		// file governs all of them; the copies are byte-identical and MSFragger is
		// invoked with the first.
		//
		// Written BEFORE the search, not after: a failed run is exactly when the
		// parameters need inspecting. They are not a resume sentinel, so their
		// presence alone never causes a step to be skipped.
		File paramsFile = MSFraggerHeadlessRunner.writeOpenSearchParams(
				OpenSearchLayout.markerDir(resultRoot, pendingSamples.get(0)), settings);
		for (int i = 1; i < pendingSamples.size(); i++) {
			Files.copy(paramsFile.toPath(), OpenSearchLayout.markerParams(resultRoot, pendingSamples.get(i)).toPath(),
					StandardCopyOption.REPLACE_EXISTING);
		}
		cb.onMessage("Marker search parameters written to each sample's " + OpenSearchLayout.MARKER_SEARCH_DIR
				+ "/ folder: " + OpenSearchLayout.SEARCH_PARAMS_FILE);

		// One invocation for every sample: the marker database is shared, so
		// MSFragger builds its fragment index once instead of per file.
		int exit = MSFraggerHeadlessRunner.run(msfraggerJar, params.getMsfraggerMaxHeap(), paramsFile, pending, cb);
		if (exit != 0) {
			throw new IOException("FragPipe-free MSFragger marker search failed (exit=" + exit + ")");
		}

		for (int i = 0; i < pending.size(); i++) {
			File produced = MSFraggerHeadlessRunner.expectedPepXml(pending.get(i));
			File sentinel = OpenSearchLayout.markerPepXml(resultRoot, pendingSamples.get(i));
			if (!produced.exists()) {
				throw new IOException("Missing MSFragger pepXML after marker search: " + produced.getAbsolutePath());
			}
			// The pin, when one was asked for, moves first: it is not the resume
			// sentinel, so a crash between the two only costs a re-run, never a
			// falsely-skipped step.
			File producedPin = MSFraggerHeadlessRunner.expectedPin(pending.get(i));
			if (producedPin.exists()) {
				Files.copy(producedPin.toPath(),
						OpenSearchLayout.markerPin(resultRoot, pendingSamples.get(i)).toPath(),
						StandardCopyOption.REPLACE_EXISTING);
			}
			if (convert) {
				File mzml = OpenSearchSpectra.collect(pending.get(i), resultRoot);
				cb.onMessage("Spectra of " + pending.get(i).getName() + " for the rest of the workflow: " + mzml);
			}
			// Copied last, so a crash mid-search leaves the sentinel absent and
			// the step correctly reruns rather than being falsely skipped.
			Files.copy(produced.toPath(), sentinel.toPath(), StandardCopyOption.REPLACE_EXISTING);
		}
	}

	/**
	 * Builds the refined FASTA for one sample from stage-1 evidence.
	 *
	 * @return {@code true} when the sample has a refined database ready for
	 *         stage-2 identification; {@code false} when stage-1 evidence is
	 *         empty and identification should be skipped for this sample
	 *         (quantification can still include it via FlashLFQ MBR)
	 */
	private static boolean extractUnmodified(DDAMagOpenSearchParameters params, File resultRoot, String[] experiment,
			MagDbItem mag, File host, boolean entrapmentRun, ProgressCallback cb) throws Exception {
		String sample = baseName(experiment[0]);
		File refined = OpenSearchLayout.refinedFasta(resultRoot, sample);
		// Staleness, not mere existence. The refined FASTA is derived from the
		// marker pepXML and - when MSBooster is on - from the edited pin, and
		// either can be rewritten by a later run without this step being asked to
		// re-run. Switching stage-1 rescoring on after a completed run is exactly
		// that case: MSBooster writes new edited pins, this step is skipped
		// because the old FASTA still exists, and stage 2 then searches a database
		// built by the very cut the rescoring was meant to replace - with nothing
		// in the output saying so.
		// The refined pepXML when Crystal-C ran, the raw one otherwise. Both the
		// staleness gate and the read below use the same file, so re-running the
		// refinement invalidates the database built from the previous masses -
		// which is the whole point: stage 1 keeps only the zero family, and an
		// uncorrected precursor puts a perfectly good PSM outside it.
		boolean corrected = round1UsesCrystalC(params);
		File markerPepXml = corrected
				? OpenSearchLayout.markerCorrectedPepXml(resultRoot, sample)
				: OpenSearchLayout.markerPepXml(resultRoot, sample);
		if (corrected && !markerPepXml.isFile()) {
			throw new IOException("Crystal-C is enabled but " + markerPepXml.getAbsolutePath() + " does not exist. "
					+ "It is written by the refinement step that runs straight after the marker searches; if that "
					+ "step was skipped, delete " + OpenSearchLayout.markerPepXml(resultRoot, sample).getName()
					+ " and re-run stage 1.");
		}
		File editedPin = OpenSearchLayout.markerEditedPin(resultRoot, sample);
		// The decoy method is checked too: a database whose decoys were made the other way is rebuilt
		// (and with it every stage-2 search, which is older than the new database).
		TargetDecoyFasta.DecoyMethod decoyMethod = params.sampleDecoyMethod();
		boolean decoysCurrent = refined.isFile() && TargetDecoyFasta.decoysMadeBy(refined, decoyMethod);
		if (OpenSearchLayout.isSatisfied(refined, markerPepXml, editedPin) && decoysCurrent) {
			// Genome selection stands; the host setting may still have changed since the database was
			// built. The genome-only targets beside it rebuild the database without repeating selection.
			// Not for entrapment: its database is the paired one OpenSearchEntrapment keeps current.
			File genomeTargets = refinedTargets(resultRoot, sample);
			if (!entrapmentRun && genomeTargets.isFile() && TargetDecoyFasta
					.ensure(OpenSearchHostDatabase.searchTargets(genomeTargets, host), refined, decoyMethod)) {
				cb.onMessage(sample + ": refined database rebuilt from the selected genomes "
						+ (host != null ? "with" : "without") + " the host proteins");
			}
			clearIdentificationSkipped(resultRoot, sample);
			cb.onMessage("Refined database already exists, skipping genome discovery: " + sample);
			return true;
		}
		if (OpenSearchLayout.isSatisfied(refined) && !decoysCurrent) {
			cb.onMessage(sample + ": refined database decoys are not " + decoyMethod.label + "; rebuilding it");
		} else if (OpenSearchLayout.isSatisfied(refined)) {
			cb.onMessage(sample + ": refined database is older than "
					+ OpenSearchLayout.newestInput(markerPepXml, editedPin).getName()
					+ "; rebuilding it so genome selection uses the current stage-1 evidence");
		}

		File unmodifiedTsv = OpenSearchLayout.unmodifiedPsms(resultRoot, sample);
		UnmodifiedPsmSelector.Result selection;

		File pepXml = markerPepXml;
		long pepXmlMb = pepXml.length() / (1024 * 1024);
		cb.onMessage(sample + ": reading marker pepXML (" + pepXmlMb + " MB) ...");
		LOGGER.info("{}: reading marker pepXML {} ({} MB)", sample, pepXml.getAbsolutePath(), pepXmlMb);

		// Rank 1 only. Stage 1 maps genomes from unmodified peptides; it does not
		// arbitrate between competing candidates, so the lower-ranked hits are
		// weight without purpose here — they inflate the read, the histogram and
		// the FDR set several-fold. Stage 2, which does arbitrate, still reads all
		// of them.
		OpenSearchPepxmlReader reader = new OpenSearchPepxmlReader(TargetDecoyFasta.DECOY_TAG);
		reader.setMaxHitRank(1);
		final ProgressCallback readCb = cb;
		List<OpenSearchPsm> psms = reader.read(pepXml, readCb::onMessage);
		reader.reportScoreNames(readCb::onMessage);
		cb.onMessage(sample + ": " + psms.size() + " top-rank PSMs read");

		cb.onMessage(sample + ": locating the zero peak (+/-" + params.toFitterConfig().zeroWindow + " Da) ...");
		selection = UnmodifiedPsmSelector.select(psms, params.toFitterConfig(), UnmodifiedPsmSelector.DEFAULT_FDR,
				readCb::onMessage);
		if (selection.unmodified.isEmpty()) {
			return markIdentificationSkipped(resultRoot, sample,
					"Marker search produced no unmodified PSMs at q<=" + UnmodifiedPsmSelector.DEFAULT_FDR
							+ "; genome discovery cannot proceed for this sample",
					cb);
		}
		// Second stage-1 pass: admit PSMs on expected mass shifts. A peptide
		// carrying a shift maps to the same genome as its unmodified form, so the
		// modification state is irrelevant to genome selection - but only shifts
		// that pass ExpectedPeakSelector's four criteria in THIS run are admitted,
		// so this stays a test of UNIMOD masses, not discovery. Open marker search
		// only: in the closed one there is no shifted peak to test (see the class).
		ExpectedPeakSelector.Config expectedConfig = params.toExpectedPeakConfig();
		ExpectedPeakSelector.Result expected = ExpectedPeakSelector.select(psms, selection, expectedConfig,
				readCb::onMessage);
		File expectedPeaksFile = OpenSearchLayout.expectedPeaks(resultRoot, sample);
		if (expectedConfig.enabled) {
			writeExpectedPeakVerdicts(expected, expectedPeaksFile);
		} else {
			// A file left by an earlier open-mode run would describe a step that did not run.
			Files.deleteIfExists(expectedPeaksFile.toPath());
		}
		cb.onMessage(String.format(
				"%s: zero-family PSM FDR (1:1 target-decoy, decoys/targets): before=%d targets/%d decoys; "
						+ "accepted at q<=%.3f=%d targets/%d decoys; estimated FDR=%.6f; calibration offset=%.5f Da",
				sample, selection.targetsBeforeFdr, selection.decoysBeforeFdr, selection.fdrThreshold,
				selection.targetsAfterFdr, selection.decoysAfterFdr, selection.estimatedFdr,
				selection.calibrationOffset));

		// --- genome selection, via the shared MAG service -------------------
		// Targets only: decoys carry no genome evidence and the selector would
		// otherwise have to strip them itself, as the FragPipe path does.
		// Stage-1 rescoring. The candidate pool is the zero FAMILY (before the
		// engine-score cut) plus the admitted expected shifts, because rescoring
		// exists precisely to re-decide where that cut falls. Null here means
		// rescoring is switched off, never that it failed: a failure throws.
		List<OpenSearchPsm> candidates = new ArrayList<>(selection.zeroFamily);
		candidates.addAll(expected.admitted);
		Stage1Rescorer.Result rescored = runStage1Rescoring(params, resultRoot, sample, candidates, cb);

		writeUnmodifiedPsms(selection, expected, rescored, unmodifiedTsv);

		List<MetaPSM> forSelection = new ArrayList<>();
		List<OpenSearchPsm> stageOne;
		if (rescored != null) {
			stageOne = new ArrayList<>(rescored.accepted);
			cb.onMessage(String.format(
					"%s: rescored stage-1 set: %d targets/%d decoys accepted at q<=%.3f (was %d targets from the "
							+ "engine-score cut)",
					sample, rescored.targetsAfterFdr, rescored.decoysAfterFdr,
					UnmodifiedPsmSelector.DEFAULT_FDR, selection.targetsAfterFdr));
		} else {
			// Reached only when stage-1 rescoring is switched off. It is not a
			// fallback: every failure inside runStage1Rescoring throws, so this
			// branch cannot be entered by anything going wrong.
			stageOne = new ArrayList<>(selection.unmodified);
			stageOne.addAll(expected.admitted);
		}
		for (OpenSearchPsm psm : stageOne) {
			if (psm.isTarget()) {
				forSelection.add(new OpenSearchMetaPsm(psm));
			}
		}
		if (forSelection.isEmpty()) {
			String detail = rescored != null
					? ("Stage-1 rescoring accepted 0 targets at q<=" + UnmodifiedPsmSelector.DEFAULT_FDR
							+ " (engine-score cut had " + selection.targetsAfterFdr + " targets)")
					: ("No target PSMs in the unmodified set after the engine-score cut (accepted "
							+ selection.targetsAfterFdr + " targets / " + selection.decoysAfterFdr + " decoys)");
			return markIdentificationSkipped(resultRoot, sample, detail, cb);
		}
		cb.onMessage(sample + ": selecting genomes from " + forSelection.size() + " unmodified target PSMs ("
				+ genomeSelectionSummary(params) + ")");

		File markerDir = OpenSearchLayout.markerDir(resultRoot, sample);
		String prefix = MagDbItem.MARKER4.equalsIgnoreCase(params.getGenomeSelectorMode()) ? "marker_genome"
				: "genome";

		MagGenomeSelectionService.Request request = new MagGenomeSelectionService.Request();
		request.magDb = mag;
		copyGenomeSelection(request, params);
		request.psms = forSelection;
		request.evidenceFile = new File(markerDir, prefix + "_evidence.tsv");
		request.rankedFile = new File(markerDir, prefix + "_ranked.tsv");
		request.selectedFile = new File(markerDir, prefix + "_selected.tsv");
		// The PSM file this step has just rewritten. Without it the service's own
		// resume check would hand back the previous run's genomes and the rescored
		// evidence would never reach the refined database.
		request.inputs = new File[] { unmodifiedTsv };
		request.maximumSelected = MagGenomeSelectionService.DEFAULT_MAXIMUM_SELECTED;
		request.inferFromProteinNamesWhenEmpty = true;
		request.allowMarker4LibraryFallback = true;
		request.logDiagnostics = MagDbItem.MARKER4.equalsIgnoreCase(params.getGenomeSelectorMode());
		request.sourceDescription = unmodifiedTsv.getAbsolutePath();

		HashSet<String> genomes = MagGenomeSelectionService.selectGenomes(request, cb);
		if (genomes.isEmpty()) {
			return markIdentificationSkipped(resultRoot, sample,
					"Genome selection returned no genomes from " + forSelection.size() + " unmodified target PSMs",
					cb);
		}
		cb.onMessage(sample + ": " + genomes.size() + " genomes selected from unmodified PSMs");

		buildRefinedFasta(mag, genomes, sample, markerDir, refined, params.sampleDecoyMethod(), host, cb);
		clearIdentificationSkipped(resultRoot, sample);
		return true;
	}

	/**
	 * After stage-1 genome selection, drop samples that have no refined database
	 * when the caller opts in; otherwise stop with the same fail-closed default
	 * DIA Deep per-file uses.
	 */
	static List<String[]> applyIdentificationGate(List<String[]> successfulExperiments,
			List<String> skippedSampleNames, PartialIdentificationHandler partialIdentificationHandler,
			ProgressCallback cb) throws Exception {
		if (skippedSampleNames == null || skippedSampleNames.isEmpty()) {
			return successfulExperiments;
		}
		String missingList = String.join(", ", skippedSampleNames);
		String detail = "Stage-1 identification skipped (no refined database) for: " + missingList;
		LOGGER.warn(detail);
		cb.onMessage(detail);

		if (successfulExperiments == null || successfulExperiments.isEmpty()) {
			throw new IOException("No samples have stage-1 genome selection results after marker search. "
					+ detail + ". Replace these files and re-run.");
		}

		if (partialIdentificationHandler != null) {
			cb.onMessage("Waiting for user decision on incomplete stage-1 results ("
					+ skippedSampleNames.size() + " skipped, " + successfulExperiments.size() + " successful)");
			boolean continueRun = partialIdentificationHandler.continueWithSuccessfulFiles(skippedSampleNames,
					successfulExperiments.size());
			if (!continueRun) {
				throw new IOException("Stopped after incomplete stage-1 results. " + detail
						+ ". Replace these files and re-run, or choose continue to identify only the "
						+ successfulExperiments.size()
						+ " successful sample(s) and still quantify every raw file together.");
			}
			cb.onMessage("Continuing identification with " + successfulExperiments.size()
					+ " sample(s); skipped samples stay in FlashLFQ quantification only: " + missingList);
			return successfulExperiments;
		}

		if (continueWithSuccessfulRequested()) {
			cb.onMessage("Continuing identification with " + successfulExperiments.size()
					+ " sample(s) (" + CONTINUE_WITH_SUCCESSFUL_PROPERTY
					+ "); no stage-2 search for: " + missingList
					+ "; all samples stay in FlashLFQ quantification");
			return successfulExperiments;
		}

		throw new IOException("Incomplete stage-1 identification results. " + detail
				+ ". Successful samples: " + successfulExperiments.size()
				+ ". Replace the missing files and re-run, or run the CLI with --continue-with-successful-files"
				+ " (GUI can offer continue-with-successful). Skipped samples then omit stage-2 identification"
				+ " but still enter FlashLFQ quantification with the rest of the cohort.");
	}

	private static boolean markIdentificationSkipped(File resultRoot, String sample, String reason,
			ProgressCallback cb) throws IOException {
		File markerDir = OpenSearchLayout.markerDir(resultRoot, sample);
		markerDir.mkdirs();
		File skipFile = OpenSearchLayout.identificationSkipped(resultRoot, sample);
		Files.writeString(skipFile.toPath(), reason == null ? "" : reason);
		LOGGER.warn("{}: skipping stage-2 identification - {}", sample, reason);
		cb.onMessage(sample + ": skipping stage-2 identification - " + reason
				+ " (written " + skipFile.getName() + ")");
		return false;
	}

	private static void clearIdentificationSkipped(File resultRoot, String sample) {
		File skipFile = OpenSearchLayout.identificationSkipped(resultRoot, sample);
		if (skipFile.isFile() && !skipFile.delete()) {
			LOGGER.warn("{}: could not delete stale {}", sample, skipFile.getAbsolutePath());
		}
	}

	/**
	 * Runs MSBooster once over every sample's marker pin, before any sample's
	 * genome selection.
	 *
	 * <h2>Why one batch rather than a call per sample</h2>
	 * DIA-NN loads its prediction model on every invocation, and
	 * {@code MagFragPipeHeadlessRunner} already documents that concurrent
	 * FragPipe instances share one set of MSBooster model files, which upstream
	 * does not support. The marker search is itself a single batched MSFragger
	 * run, so there is a natural point where every pin exists at once: MSBooster
	 * is invoked here, with all of them, exactly as FragPipe invokes it once over
	 * every run in its manifest.
	 *
	 * <h2>File layout</h2>
	 * Modelled on FragPipe's: the parameter file and the shared {@code MSBooster/}
	 * prediction folder sit at the search-results root, and each
	 * {@code <sample>_edited.pin} lands beside the {@code <sample>.pin} it was
	 * derived from. Nothing is renamed - MSBooster pairs a pin with its spectra
	 * file by base name, and a suffixed pin fails with an error that blames the
	 * mzML directory instead.
	 *
	 * <h2>No silent degradation</h2>
	 * If MSBooster is switched on, every part of it must work: a missing pin, a
	 * missing jar or DIA-NN, a non-zero exit, an absent edited pin or one without
	 * the expected feature columns all abort the run. A stage that quietly
	 * produces a lesser result is indistinguishable in its output from one that
	 * worked - which is how the 2026-09-10 12mix benchmark completed all nine
	 * runs on the engine-score cut without anyone noticing. Deliberate opt-out is
	 * {@code stage1_msbooster_enabled=false}, not failure.
	 */
	private static void runStage1MsBooster(DDAMagOpenSearchParameters params, File resultRoot,
			List<String[]> experiments, ProgressCallback cb) throws IOException {
		if (!params.isStage1RescoringEnabled() || !params.isStage1MsboosterEnabled()) {
			return;
		}
		File jar = params.getMsboosterJar();
		File diann = params.getDiannExecutable();
		if (jar == null || !jar.isFile()) {
			throw new IOException("Stage-1 MSBooster is enabled but MSBooster.jar was not found. Set "
					+ "Resources > FragPipe to an installation whose tools folder contains it.");
		}
		if (diann == null || !diann.isFile()) {
			throw new IOException("Stage-1 MSBooster is enabled but DIA-NN was not found. MSBooster predicts "
					+ "through DIA-NN and needs the 1.8.x prediction interface that FragPipe bundles under "
					+ "tools/diann/; a standalone 2.x build will not work.");
		}

		List<File> pins = new ArrayList<>();
		List<String> pending = new ArrayList<>();
		File fraggerParams = null;
		for (String[] experiment : experiments) {
			String sample = baseName(experiment[0]);
			// Resume first, then inputs. The edited pin is this step's own
			// sentinel; once it exists the raw .pin has done its job and nothing
			// downstream reads it again, so a sample that is already finished must
			// not fail because its input was tidied away - these pins are ~54 MB
			// each and are exactly the kind of file a disk cleanup removes.
			// Older than the marker pepXML means the search re-ran (for example
			// after a change of round-1 settings) and the pin predates it.
			if (OpenSearchLayout.isSatisfied(OpenSearchLayout.markerEditedPin(resultRoot, sample),
					OpenSearchLayout.markerPepXml(resultRoot, sample))) {
				cb.onMessage("MSBooster output already exists, skipping: " + sample);
				continue;
			}
			File pin = OpenSearchLayout.markerPin(resultRoot, sample);
			if (!pin.isFile()) {
				throw new IOException("Stage-1 MSBooster is enabled but there is no .pin for " + sample + ": "
						+ pin.getAbsolutePath() + ", and no " + OpenSearchLayout.markerEditedPin(resultRoot, sample)
								.getName()
						+ " either. Round one must run with output_format = pepXML_pin; if the pin was deleted "
						+ "after a completed run, delete this sample's marker pepXML too so the search re-runs.");
			}
			pins.add(pin);
			pending.add(sample);
			if (fraggerParams == null) {
				fraggerParams = OpenSearchLayout.markerParams(resultRoot, sample);
			}
		}
		if (pins.isEmpty()) {
			return;
		}
		cb.onMessage("Stage-1 MSBooster: predicting spectra and retention times for " + pins.size()
				+ " sample(s) in one invocation: " + String.join(", ", pending));

		MSBoosterHeadlessRunner.Settings settings = new MSBoosterHeadlessRunner.Settings();
		settings.diann = diann;
		settings.unimodObo = params.getUnimodObo();
		settings.threads = params.getThreadCount();
		// Round one's parameters, so MSBooster reads the same tolerances and
		// modifications the search used instead of falling back to its defaults.
		// Every sample's copy is byte-identical; the first is as good as any.
		settings.fraggerParams = fraggerParams;
		settings.pinFiles = pins;
		settings.mzmlDirectories = MSBoosterHeadlessRunner.listOf(readableSpectraDir(params.getExperimentData(),
				resultRoot));

		File msboosterParams = MSBoosterHeadlessRunner.writeParams(
				OpenSearchLayout.searchResultsRoot(resultRoot), settings);
		int exit;
		try {
			exit = MSBoosterHeadlessRunner.run(jar, params.getMsboosterMaxHeap(), msboosterParams, cb);
		} catch (InterruptedException e) {
			Thread.currentThread().interrupt();
			throw new IOException("Interrupted while running MSBooster", e);
		}
		if (exit != 0) {
			throw new IOException("MSBooster exited " + exit + ". Its parameters are in "
					+ msboosterParams.getAbsolutePath() + "; the DIA-NN it was pointed at is "
					+ diann.getAbsolutePath() + " (MSBooster needs the 1.8.x prediction interface, not 2.x).");
		}
		for (String sample : pending) {
			File edited = OpenSearchLayout.markerEditedPin(resultRoot, sample);
			if (!edited.isFile()) {
				throw new IOException("MSBooster reported success but wrote no " + edited.getName() + " for "
						+ sample + ". Check its output above and the parameters in "
						+ msboosterParams.getAbsolutePath() + ".");
			}
		}
		// MSBooster names its own output folder and always calls it "MSBooster",
		// so a later stage-2 batch over sample_db pins would resolve to the same
		// place and overwrite these predictions. Move them somewhere stage-named
		// before that can happen.
		File predictions = MSBoosterHeadlessRunner.relocatePredictions(pins,
				OpenSearchLayout.msboosterDir(resultRoot, OpenSearchLayout.MARKER_SEARCH_DIR));
		MSBoosterHeadlessRunner.writeRunNote(predictions, pending);
		if (predictions != null) {
			File plots = new File(predictions, MSBoosterHeadlessRunner.PLOTS_DIR);
			for (String sample : pending) {
				MSBoosterHeadlessRunner.copySamplePlots(plots, sample,
						OpenSearchLayout.markerDir(resultRoot, sample));
			}
		}
		cb.onMessage("Stage-1 MSBooster: edited pins and QC plots written beside each marker pin; the shared "
				+ "prediction library is in "
				+ (predictions == null ? "the MSBooster folder" : predictions.getAbsolutePath()));
	}

	/**
	 * Rescores this sample's stage-1 candidates, reading MSBooster's predicted
	 * features from the edited pin {@link #runStage1MsBooster} produced.
	 * <p>
	 * <b>No silent degradation.</b> With MSBooster switched on, a missing edited
	 * pin, a pin without the expected feature columns, or a classifier failure
	 * all abort the run rather than falling back to the engine-score cut. See
	 * {@link #runStage1MsBooster} for why.
	 */
	private static Stage1Rescorer.Result runStage1Rescoring(DDAMagOpenSearchParameters params, File resultRoot,
			String sample, List<OpenSearchPsm> candidates, ProgressCallback cb) throws IOException {
		if (!params.isStage1RescoringEnabled()) {
			return null;
		}
		Map<String, double[]> predicted = null;
		List<String> predictedNames = null;
		if (params.isStage1MsboosterEnabled()) {
			File edited = OpenSearchLayout.markerEditedPin(resultRoot, sample);
			if (!edited.isFile()) {
				throw new IOException("Stage-1 MSBooster is enabled but " + edited.getAbsolutePath()
						+ " does not exist. It is written by the batched MSBooster step that runs after the "
						+ "marker searches; if that step was skipped, the pin beside it is stale - delete "
						+ OpenSearchLayout.markerPin(resultRoot, sample).getName() + " and re-run stage 1.");
			}
			MSBoosterPin.Features features = MSBoosterPin.readFeatures(edited,
					MSBoosterHeadlessRunner.DEFAULT_FEATURE_COLUMNS);
			if (features.isEmpty()) {
				throw new IOException("MSBooster produced " + edited.getName() + " but it carries none of the "
						+ "expected feature columns " + MSBoosterHeadlessRunner.DEFAULT_FEATURE_COLUMNS
						+ ". Missing: " + features.missing);
			}
			predicted = features.byKey;
			predictedNames = features.names;
			// Only the candidates this stage consumes are looked up, so the rows
			// the edited pin holds for the rest of the pin cost nothing beyond
			// this read. Reported all the same: a collapse in the matched count
			// between runs is the first sign the pin and pepXML have drifted
			// apart.
			// Every candidate must be present. Round one runs at
			// output_report_topN = 1, so the pepXML and the pin describe the same
			// rank-1 PSM for the same spectrum and the correspondence is exact. A
			// candidate with no row would be scored on search features while its
			// neighbours were scored on search plus prediction - not a smaller
			// result but an incoherent one, and invisible in the output. If this
			// ever fires legitimately the cause is worth knowing, not absorbing.
			int matched = 0;
			OpenSearchPsm firstMissing = null;
			for (OpenSearchPsm psm : candidates) {
				if (predicted.containsKey(MSBoosterPin.key(psm.getFileName(), psm.getScan(), psm.getHitRank()))) {
					matched++;
				} else if (firstMissing == null) {
					firstMissing = psm;
				}
			}
			if (matched < candidates.size()) {
				throw new IOException(String.format(
						"%d of %d stage-1 candidates have no row in %s (first: %s scan %d rank %d). The pepXML and "
								+ "the pin disagree about which PSMs the search reported; they must come from the "
								+ "same round-one search, written with output_format = pepXML_pin.",
						candidates.size() - matched, candidates.size(), edited.getName(),
						firstMissing.getFileName(), firstMissing.getScan(), firstMissing.getHitRank()));
			}
			cb.onMessage(String.format("%s: MSBooster supplied %d predicted features for all %d stage-1 candidates",
					sample, features.names.size(), candidates.size()));
		}
		return Stage1Rescorer.rescore(candidates, params.toStage1RescorerConfig(), predicted, predictedNames,
				cb::onMessage);
	}

	/**
	 * Writes the selected genomes' sequences to a targets-only FASTA, then hands
	 * it to {@link TargetDecoyFasta} for the 1:1 {@code rev_} decoy pass — rather
	 * than interleaving decoys while writing, as the FragPipe path does. Reusing
	 * the shared builder means the refined database gets the same
	 * one-decoy-per-target validation the marker database does, which is what the
	 * downstream {@code decoys/targets} FDR estimator assumes.
	 */
	private static void buildRefinedFasta(MagDbItem mag, HashSet<String> genomes, String sample, File markerDir,
			File refined, TargetDecoyFasta.DecoyMethod decoyMethod, File host, ProgressCallback cb) throws IOException {
		File originalDbFolder = mag.getOriginalDbFolder();
		if (originalDbFolder == null || !originalDbFolder.exists()) {
			throw new IOException("Original DB folder not found in MAG DB: " + originalDbFolder);
		}

		File targetsOnly = new File(markerDir, sample + REFINED_TARGETS_SUFFIX);
		cb.onMessage(sample + ": writing refined database from " + genomes.size() + " genomes ...");
		int written = 0;
		int genomesDone = 0;
		try (PrintWriter writer = new PrintWriter(targetsOnly, "UTF-8")) {
			for (String genome : genomes) {
				File faa = new File(originalDbFolder, genome + ".faa");
				File fa = new File(originalDbFolder, genome + ".fasta");
				File src = faa.exists() ? faa : (fa.exists() ? fa : null);
				if (src == null) {
					continue;
				}
				for (String line : Files.readAllLines(src.toPath())) {
					writer.println(line);
					if (line.startsWith(">")) {
						written++;
					}
				}
				if (++genomesDone % 100 == 0) {
					cb.onMessage(sample + ": " + genomesDone + "/" + genomes.size() + " genomes written ("
							+ written + " sequences)");
				}
			}
		}
		if (written == 0) {
			throw new IOException("Refined database for " + sample + " is empty: no sequences found for "
					+ genomes.size() + " selected genomes under " + originalDbFolder.getAbsolutePath());
		}

		File searched = OpenSearchHostDatabase.searchTargets(targetsOnly, host);
		TargetDecoyFasta.ensure(searched, refined, decoyMethod);
		int hostProteins = searched == targetsOnly ? 0 : TargetDecoyFasta.count(refined).targets - written;
		cb.onMessage(sample + ": refined database built (" + written + " genome target sequences"
				+ (host != null ? " + " + hostProteins + " host" : "") + " + 1:1 " + decoyMethod.label + " decoys)");
	}

	static final String REFINED_TARGETS_SUFFIX = "_refined_targets.fasta";

	/** The selected genomes' proteins, targets only, that a sample's refined database is built from. */
	static File refinedTargets(File resultRoot, String sample) {
		return new File(OpenSearchLayout.markerDir(resultRoot, sample), sample + REFINED_TARGETS_SUFFIX);
	}

	private static List<File> runSampleDbSearches(DDAMagOpenSearchParameters params, File msfraggerJar, File resultRoot,
			List<String[]> experiments, ProgressCallback cb) throws Exception {
		List<File> pepxmls = new ArrayList<>();
		for (String[] experiment : experiments) {
			String sample = baseName(experiment[0]);
			File sentinel = OpenSearchLayout.sampleDbPepXml(resultRoot, sample);
			// Same staleness rule as step B, one link further down the chain: a
			// rebuilt refined FASTA must invalidate the stage-2 search that was
			// run against the old one, or the upgrade stops at the database and
			// never reaches the results.
			File refined = OpenSearchLayout.refinedFasta(resultRoot, sample);

			// Staleness, by configuration rather than by timestamp - mirrors the
			// Sage round-2 search below. A changed output_report_topN (or anything
			// else baseSettings derives from params) leaves the previous pepXML
			// intact and newer than the refined FASTA, so timestamp comparison
			// alone silently keeps serving a search that ran under the old
			// settings. Found 2026-09-17: exactly this left a topN=1 pepXML on
			// disk after the rank policy was reset to 5, and a run three days
			// later silently reused it with no message saying so. Only checked
			// when the refined FASTA exists; when it does not, the explicit check
			// below reports that clearly instead.
			MSFraggerHeadlessRunner.OpenSearchSettings settings = null;
			String wantedParams = null;
			File paramsFile = OpenSearchLayout.sampleDbParams(resultRoot, sample);
			if (refined.exists()) {
				settings = baseSettings(params, refined);
				// Round two reports one candidate per spectrum by default. The
				// extra ranks were justified as peak-discovery votes and
				// classifier negatives, but measured on 12mix they buy access to
				// a population at roughly 50-60% genome-attribution purity: at
				// every matched entrapment budget, rank-1 alone recovered more
				// true PSMs. Negatives are not scarce either - 40.9% of
				// candidates outside any peak are decoys.
				settings.outputReportTopN = Math.max(1, params.getRound2OutputReportTopN());
				wantedParams = MSFraggerHeadlessRunner.previewOpenSearchParamsText(settings);
			}
			boolean paramsUnchanged = wantedParams != null && paramsFile.exists()
					&& wantedParams.trim().equals(readFileText(paramsFile).trim());

			if (paramsUnchanged && OpenSearchLayout.isSatisfied(sentinel, refined)) {
				cb.onMessage("Sample-specific search already exists under the same configuration, skipping: "
						+ sample);
				pepxmls.add(sentinel);
				continue;
			}
			if (OpenSearchLayout.isSatisfied(sentinel)) {
				cb.onMessage(sample + ": " + (wantedParams != null && !paramsUnchanged
						? "the MSFragger configuration has changed since this search ran"
						: "sample-specific search is older than its refined database")
						+ "; re-running it");
			}

			if (!refined.exists()) {
				throw new IOException("Refined database missing for " + sample + ": " + refined.getAbsolutePath()
						+ ". Genome selection must run before the stage-2 search.");
			}
			File sampleDbDir = OpenSearchLayout.sampleDbDir(resultRoot, sample);
			sampleDbDir.mkdirs();

			MSFraggerHeadlessRunner.writeOpenSearchParams(sampleDbDir, settings);

			// Per-sample invocation: each sample has its own refined database, so
			// the fragment-index reuse that stage 1 exploits is unavailable.
			int exit = MSFraggerHeadlessRunner.run(msfraggerJar, params.getMsfraggerMaxHeap(), paramsFile,
					Collections.singletonList(new File(experiment[0])), cb);
			if (exit != 0) {
				throw new IOException("MSFragger sample-specific search failed for " + sample + " (exit=" + exit + ")");
			}
			File produced = MSFraggerHeadlessRunner.expectedPepXml(new File(experiment[0]));
			if (!produced.exists()) {
				throw new IOException("Missing MSFragger pepXML after sample-specific search: " + produced.getAbsolutePath());
			}
			Files.copy(produced.toPath(), sentinel.toPath(), StandardCopyOption.REPLACE_EXISTING);
			pepxmls.add(sentinel);
		}
		return pepxmls;
	}

	/**
	 * Closed stage of the closed-then-open identification: each sample's refined database searched
	 * closed ({@link #closedStageSettings}) into {@code sample_db_closed/}. Empty when the stage is
	 * off. Staleness as for the open search: by parameter text and against the refined FASTA.
	 */
	private static List<File> runSampleDbClosedSearches(DDAMagOpenSearchParameters params, File msfraggerJar,
			File resultRoot, List<String[]> experiments, ProgressCallback cb) throws Exception {
		List<File> pepxmls = new ArrayList<>();
		if (!params.isClosedStageEnabled()) {
			return pepxmls;
		}
		for (String[] experiment : experiments) {
			String sample = baseName(experiment[0]);
			File sentinel = OpenSearchLayout.sampleDbClosedPepXml(resultRoot, sample);
			File refined = OpenSearchLayout.refinedFasta(resultRoot, sample);
			if (!refined.exists()) {
				throw new IOException("Refined database missing for " + sample + ": " + refined.getAbsolutePath()
						+ ". Genome selection must run before the closed search.");
			}
			MSFraggerHeadlessRunner.OpenSearchSettings settings = closedStageSettings(params, refined);
			String wantedParams = MSFraggerHeadlessRunner.previewOpenSearchParamsText(settings);
			File paramsFile = OpenSearchLayout.sampleDbClosedParams(resultRoot, sample);
			boolean paramsUnchanged = paramsFile.exists()
					&& wantedParams.trim().equals(readFileText(paramsFile).trim());
			if (paramsUnchanged && OpenSearchLayout.isSatisfied(sentinel, refined)) {
				cb.onMessage("Closed search already exists under the same configuration, skipping: " + sample);
				pepxmls.add(sentinel);
				continue;
			}
			if (OpenSearchLayout.isSatisfied(sentinel)) {
				cb.onMessage(sample + ": " + (!paramsUnchanged
						? "the closed-search configuration has changed since this search ran"
						: "closed search is older than its refined database")
						+ "; re-running it");
			}
			File closedDir = OpenSearchLayout.sampleDbClosedDir(resultRoot, sample);
			closedDir.mkdirs();
			MSFraggerHeadlessRunner.writeOpenSearchParams(closedDir, settings);
			cb.onMessage("Step C: closed search of the refined database - " + sample);
			int exit = MSFraggerHeadlessRunner.run(msfraggerJar, params.getMsfraggerMaxHeap(), paramsFile,
					Collections.singletonList(new File(experiment[0])), cb);
			if (exit != 0) {
				throw new IOException("MSFragger closed search failed for " + sample + " (exit=" + exit + ")");
			}
			File produced = MSFraggerHeadlessRunner.expectedPepXml(new File(experiment[0]));
			if (!produced.exists()) {
				throw new IOException("Missing MSFragger pepXML after closed search: " + produced.getAbsolutePath());
			}
			Files.copy(produced.toPath(), sentinel.toPath(), StandardCopyOption.REPLACE_EXISTING);
			pepxmls.add(sentinel);
		}
		return pepxmls;
	}

	/**
	 * @param closedPsms the closed search's PSMs for the closed-then-open identification, or null
	 *                   for the open-only workflow
	 */
	/** Writes {@link OpenSearchLayout#identificationMode} when its content changes, and returns it. */
	static File writeIdentificationModeStamp(DDAMagOpenSearchParameters params, File resultRoot, String engineMode)
			throws IOException {
		return writeIdentificationModeStamp(params, resultRoot, engineMode, null);
	}

	/**
	 * @param spectraFiles the runs' spectrum files, for recording which have an mzML the
	 *                     predicted-spectrum features can read; null records none
	 */
	static File writeIdentificationModeStamp(DDAMagOpenSearchParameters params, File resultRoot, String engineMode,
			List<File> spectraFiles) throws IOException {
		File stamp = OpenSearchLayout.identificationMode(resultRoot, engineMode);
		String mode = !params.isClosedStageEnabled() ? "open_only"
				: params.isJointIdentification()
						? "joint_closed_open closed_stage_precursor_ppm=" + params.getClosedStagePrecursorPpm()
								+ " closed_stage_topn=" + params.getClosedStageTopN() + " ms2_features="
								+ params.isJointMs2Features() + ms2InputsStamp(params, spectraFiles)
								+ " ms2_max_open_rank=" + params.getJointMs2MaxOpenRank()
								+ " fragment_ppm=" + params.getJointFragmentTolerancePpm() + " accept_unexplained="
								+ params.isJointAcceptUnexplained() + " compare_arms=" + params.isJointCompareArms()
								+ " withheld=" + params.getJointWithheldFeatures()
						: "closed_then_open closed_stage_precursor_ppm=" + params.getClosedStagePrecursorPpm();
		if (params.isClosedStageEnabled() && (params.isCascadeMs2Features() || params.isCascadeClosedRerank())) {
			// Only when on, so runs made before these options existed stay current.
			mode += " cascade_closed_rerank=" + params.isCascadeClosedRerank() + " cascade_ms2_features="
					+ params.isCascadeMs2Features()
					+ (params.isCascadeMs2Features() ? cascadeMs2Stamp(params, spectraFiles) : "");
		}
		if (params.isClosedStageEnabled() && params.isJointIdentification() && params.isIdentificationV2()) {
			mode += " identification_version=v2 model=" + identificationV2Model().describe()
					+ " shift_localisation=" + new ShiftLocaliser(params.getJointFragmentTolerancePpm());
		}
		if (!stamp.isFile() || !mode.equals(readFileText(stamp).trim())) {
			stamp.getParentFile().mkdirs();
			Files.write(stamp.toPath(), (mode + System.lineSeparator()).getBytes(java.nio.charset.StandardCharsets.UTF_8));
		}
		return stamp;
	}

	/**
	 * What the predicted-spectrum features can actually use: whether the DIA-NN executable exists
	 * and how many runs have an mzML. The features switch themselves off without them, so a run
	 * made before DIA-NN was installed would otherwise count as current after it is.
	 */
	private static String ms2InputsStamp(DDAMagOpenSearchParameters params, List<File> spectraFiles) {
		if (!params.isJointMs2Features()) {
			return "";
		}
		return cascadeMs2Stamp(params, spectraFiles);
	}

	private static String cascadeMs2Stamp(DDAMagOpenSearchParameters params, List<File> spectraFiles) {
		File diann = params.getDiannExecutable();
		boolean diannReady = diann != null && diann.isFile();
		int mzml = spectraFiles == null ? 0 : VariantFormRanking.mzmlBySample(spectraFiles).size();
		return " ms2_inputs=" + (diannReady ? "diann:" + diann.getAbsolutePath() : "no_diann") + ",mzml_runs:" + mzml;
	}

	private static ModificationValidator.Result validate(DDAMagOpenSearchParameters params, List<OpenSearchPsm> psms,
			List<OpenSearchPsm> closedPsms, List<File> spectraFiles, File resultRoot, ProgressCallback cb,
			AtomicReference<SubstitutionReassigner.Result> substitutionsOut) throws Exception {
		boolean cascade = closedPsms != null;
		boolean joint = cascade && params.isJointIdentification();
		// The cascade's closed stage and the open-only classifier are needed by the joint path only
		// to write the comparison arms.
		boolean compareArms = joint && params.isJointCompareArms();
		SemiSupervisedPsmClassifier.Config classifierConfig = new SemiSupervisedPsmClassifier.Config();
		classifierConfig.folds = params.getClassifierFolds();
		classifierConfig.maxIterations = params.getClassifierMaxIterations();
		classifierConfig.trainQValue = params.getClassifierTrainQValue();

		ClosedOpenCascade.ClosedStage closed = null;
		if (cascade && (!joint || compareArms)) {
			cb.onMessage("Closed stage: scoring " + closedPsms.size() + " closed-search PSMs with a classifier "
					+ "trained on them alone ...");
			String[] closedExtra = new String[0];
			if (params.isCascadeMs2Features()
					&& attachCascadePredictions(params, closedPsms, spectraFiles, resultRoot, "closed", cb) > 0) {
				closedExtra = PredictedSpectrumFeatures.PSM_FEATURES;
			}
			closed = ClosedOpenCascade.validateClosed(closedPsms, classifierConfig, CASCADE_FDR, cb::onMessage,
					params.isCascadeClosedRerank(), closedExtra);
			cb.onMessage(String.format(java.util.Locale.ROOT,
					"Closed stage: %,d of %,d spectra accepted at %.0f%% FDR (%,d targets, %,d decoys); the open "
							+ "search identifies only the rest",
					closed.accepted.size(), closed.winners.size(), CASCADE_FDR * 100, closed.acceptedTargets,
					closed.acceptedDecoys));
		}

		double[] deltaMasses = new double[psms.size()];
		boolean[] targetsOnly = new boolean[psms.size()];
		double[] discoveryWeights = new double[psms.size()];
		String[] fileGroups = new String[psms.size()];
		Map<String, Integer> candidatesPerSpectrumAndClass = new HashMap<>();
		for (int i = 0; i < psms.size(); i++) {
			OpenSearchPsm psm = psms.get(i);
			deltaMasses[i] = psm.getMassDiff();
			targetsOnly[i] = psm.isTarget();
			fileGroups[i] = psm.getFileName();
			// Counted over RANK-1 candidates only, matching the detection population
			// below. Weighting over all 5 ranks while detecting on rank 1 alone would
			// leave each spectrum contributing 1/topN of a vote instead of 1, and
			// minPeakArea and kdeMinFileSupport are absolute thresholds on that
			// weighted support - so every peak would silently need topN times more
			// spectra to be found.
			//
			// One vote per spectrum per target class, however many engines ran: one
			// engine gives its rank-1 row weight 1; two engines give each rank-1 row
			// 0.5, so agreement puts the whole vote on one delta mass and
			// disagreement splits it. Identical behaviour at any engine count.
			if (psm.getHitRank() == 1) {
				String key = psm.getFileName() + "\u0000" + psm.getScan() + "\u0000" + psm.isTarget();
				candidatesPerSpectrumAndClass.merge(key, 1, Integer::sum);
			}
		}
		for (int i = 0; i < psms.size(); i++) {
			OpenSearchPsm psm = psms.get(i);
			if (psm.getHitRank() != 1) {
				discoveryWeights[i] = 0.0;
				continue;
			}
			String key = psm.getFileName() + "\u0000" + psm.getScan() + "\u0000" + psm.isTarget();
			discoveryWeights[i] = 1.0 / candidatesPerSpectrumAndClass.get(key);
		}

		final ProgressCallback fitCb = cb;

		// Every file contributes to peak picking. The discovery/gate file split was
		// removed on 2026-09-13 along with the significance test it protected: a
		// threshold on a measured decoy fraction is not a hypothesis test, so there
		// is nothing to hold out, and holding files back cost half the discovery
		// data and made a modification confined to the held-out samples invisible.
		//
		// Detection uses RANK 1 ONLY. Before rescoring, the engine's own ranking is
		// the only ordering available, so admitting ranks 2-5 here would contradict
		// the evidence the step depends on. Measured on 12mix: 57.4% of rank-1
		// candidates land inside a fitted peak against 16.9% of ranks 2-5, and under
		// per-spectrum weighting those lower ranks were carrying roughly 80% of each
		// spectrum's vote into diffuse background.
		boolean[] rankOne = new boolean[psms.size()];
		int rankOneCount = 0;
		for (int i = 0; i < psms.size(); i++) {
			rankOne[i] = psms.get(i).getHitRank() == 1;
			if (rankOne[i]) {
				rankOneCount++;
			}
		}
		cb.onMessage("Peak picking on " + rankOneCount + " rank-1 candidates from all "
				+ new TreeSet<>(Arrays.asList(fileGroups)).size() + " file(s)");

		DeltaMassPeakFitter.Config fitterConfig = params.toFitterConfig();
		fitterConfig.detectionRows = rankOne;
		DeltaMassPeakFitter.Result fit = DeltaMassPeakFitter.fit(deltaMasses, targetsOnly, discoveryWeights,
				targetsOnly, fileGroups, fitterConfig, fitCb::onMessage);
		writePeaks(fit, psms, OpenSearchLayout.fittedPeaks(resultRoot, params.getSearchEngines()));

		OpenSearchFeatures.applyPeakContext(psms, fit);

		// Is there fragment-level evidence that any delta mass is actually on the
		// peptide? MSFragger scored every spectrum with and without the shift
		// (localize_delta_mass = 1), so this costs nothing to ask. Peaks with no
		// support are precursor errors, not modifications, and their PSMs are
		// folded onto the zero family as the unmodified peptides they are.
		// Runs before the filter and before scoring, so every later step sees the
		// corrected assignment.
		DeltaMassEvidenceGate.Config evidenceConfig = params.toDeltaMassEvidenceConfig();
		evidenceConfig.foldRejected = !cascade;
		DeltaMassEvidenceGate.Result fragmentEvidence = DeltaMassEvidenceGate.apply(psms, fit, evidenceConfig,
				fitCb::onMessage);

		// Is the delta explained by the SEQUENCE being wrong rather than the mass?
		// A single amino-acid substitution is invisible to target-decoy - the wrong
		// peptide is itself a target - and survives the fragment-evidence test
		// above, because a substituted residue really does shift part of the ion
		// series. Only the substitution mass table plus the localization can see it.
		TrypticPeptideIndex peptideIndex = buildPeptideIndex(params, resultRoot, cb);
		SubstitutionReassigner.Result substitutions = SubstitutionReassigner.analyse(psms,
				peptideIndex, params.toSubstitutionConfig(), fitCb::onMessage);
		if (substitutionsOut != null) {
			substitutionsOut.set(substitutions);
		}

		// The filter runs here - before any PSM is scored - on rank-1 candidates.
		// Being score-independent is what lets its verdict serve as a classifier
		// feature below without circularity.
		ModificationValidator.Gate gate = ModificationValidator.filterPeaks(psms,
				params.toValidatorConfig(), fragmentEvidence.supportedPeakIds, fitCb::onMessage);
		// Withdrawn, not rejected: these peaks may pass the screen perfectly well;
		// the point is that they are not modification hypotheses at all, so
		// counting them as modifications would be wrong however clean they are.
		// A wrong peptide with a consistent offset passes the screen and entrapment
		// alike; on 12mix the +/-14 peaks do, with 18-20% of PSMs outside the true genomes.
		// Localisation decides: where the peak's localised PSMs sit (PeakSiteShare).
		Set<Integer> screened = new HashSet<>(gate.acceptedPeakIds);
		for (ModificationValidator.PeakVerdict verdict : gate.verdicts) {
			if (verdict.massClass == ModificationValidator.MassClass.ZERO) {
				screened.remove(verdict.peakId);
			}
		}
		PeakSiteShare.Result siteShare = PeakSiteShare.measure(psms, screened,
				new ShiftLocaliser(params.getJointFragmentTolerancePpm()), VariantFormRanking.mzmlBySample(spectraFiles),
				new PeakSiteShare.Config(), fitCb::onMessage);
		PeakSiteShare.write(siteShare, OpenSearchLayout.peakSiteShares(resultRoot, params.getSearchEngines()));
		gate = gate.withoutPeaks(siteShare.withdrawnPeakIds());
		if (!substitutions.peaks.isEmpty()) {
			SubstitutionReassigner.writePeakReport(substitutions,
					OpenSearchLayout.substitutionPeaks(resultRoot, params.getSearchEngines()));
		}
		if (!substitutions.proposals.isEmpty()) {
			File candidates = OpenSearchLayout.substitutionCandidates(resultRoot, params.getSearchEngines());
			int written = SubstitutionReassigner.writeCandidateFasta(substitutions, candidates);
			cb.onMessage("Substitution candidates: " + written + " proposed sequence(s) written to "
					+ candidates.getName() + ". Searching these against the affected spectra is what decides "
					+ "whether each PSM is a modification or a misassigned sequence; until then they are "
					+ "reported as identifications, not modifications.");
		}
		for (OpenSearchPsm psm : psms) {
			psm.setInValidatedPeak(psm.getPeakId() >= 0 && gate.acceptedPeakIds.contains(psm.getPeakId()));
		}

		OpenSearchFeatures.computeAggregateCounts(psms);

		if (joint && !compareArms) {
			ModificationValidator.Result validation = runJoint(params, psms, closedPsms, gate, spectraFiles,
					resultRoot, classifierConfig, cb, null);
			alignSubstitutionVerdicts(params, validation, peptideIndex, cb);
			writePotentialModifications(fit, substitutions, validation,
					OpenSearchLayout.potentialModifications(resultRoot, params.getSearchEngines()));
			return validation;
		}

		if (cascade && params.isCascadeMs2Features()) {
			// The open classifier reads the block when it is there (CombinedPsmFeatures.PREDICTED).
			attachCascadePredictions(params, psms, spectraFiles, resultRoot, "open", cb);
		}

		if (params.isRestrictPositivesToValidatedPeaks()) {
			// Learn from PSMs that are both confident and explained. Scoring still
			// covers every candidate; only the seed set is narrowed.
			classifierConfig.positiveFilter = OpenSearchPsm::isInValidatedPeak;
		}

		// One instance per (spectrum, peptide) carrying both engines' evidence, or
		// - when only one engine ran, or the combined path is switched off - one
		// instance per engine row, calibrated onto a common scale afterwards.
		//
		// The combined path exists because that calibration is not achieving what
		// it was designed to achieve: `-log10(engine-local q)` is comparable
		// across engines only if each engine's q-value is accurate, and on the
		// 2026-09-10 run Sage's is optimistic - at equal calibrated score its rows
		// are about 1.8x more likely to be decoys. Merging first removes both the
		// engine-local q-value and the scale that had to be reconciled.
		List<OpenSearchPsm> collapseInput;
		String engineDiagnostic;
		if (params.usesCombinedPsmInstances()) {
			int[] unkeyable = new int[1];
			List<CombinedPsmCandidate> candidates = CombinedPsmCandidate.merge(psms, unkeyable);
			if (unkeyable[0] > 0) {
				cb.onMessage("Combined candidates: " + unkeyable[0]
						+ " rows had no peptide or file name and could not be merged");
			}
			CombinedPsmClassifier.Result scored = CombinedPsmClassifier.classify(candidates, classifierConfig,
					fitCb::onMessage);
			cb.onMessage(String.format(java.util.Locale.ROOT,
					"Combined candidates: %,d engine rows -> %,d candidates, %,d (%.1f%%) proposed by both engines",
					psms.size(), candidates.size(), scored.bothEngineCount,
					100.0 * scored.bothEngineCount / Math.max(1, candidates.size())));
			collapseInput = new ArrayList<>(candidates.size());
			for (int i = 0; i < candidates.size(); i++) {
				CombinedPsmCandidate candidate = candidates.get(i);
				// The score belongs to the candidate, so it is written onto every
				// engine row that produced it - otherwise scored_psms.tsv would show
				// a blank score for the row that did not happen to be the primary.
				String support = candidate.foundByBoth() ? "BOTH"
						: candidate.foundByMsfragger() ? "MSFRAGGER" : "SAGE";
				applyScore(candidate.msfragger(), scored.scores[i], scored.logits[i], support);
				applyScore(candidate.sage(), scored.scores[i], scored.logits[i], support);
				// One row per candidate goes downstream. Sending both would let a
				// peptide contest its own spectrum twice, once per engine, which is
				// double-counting rather than extra evidence.
				collapseInput.add(candidate.primary());
			}
			engineDiagnostic = "combined candidates, one model";
		} else {
			PerEngineClassifier.Result scored = PerEngineClassifier.classify(psms, classifierConfig,
					fitCb::onMessage);
			collapseInput = psms;
			StringBuilder perEngine = new StringBuilder();
			for (Map.Entry<OpenSearchPsm.SearchEngine, Integer> e : scored.acceptedByEngine.entrySet()) {
				if (perEngine.length() > 0) {
					perEngine.append(", ");
				}
				perEngine.append(e.getKey()).append(" ").append(e.getValue()).append("/")
						.append(scored.countByEngine.get(e.getKey()));
			}
			engineDiagnostic = perEngine.toString();
		}

		// Collapse to one PSM per spectrum before any FDR is computed. Target-decoy
		// FDR is defined over spectra, not candidates: with output_report_topN > 1
		// a single spectrum otherwise contributes several rows, so decoy counts
		// stop being independent, the accepted count reports candidates rather
		// than identifications, and two candidates from one spectrum could both
		// be accepted. Percolator and MS2Rescore collapse at exactly this point.
		//
		// This is not a filter on evidence: the classifier has already scored
		// every candidate and was free to promote a lower-ranked one over the
		// engine's top hit. Only the unit of analysis is being fixed.
		SpectrumLevelCollapse.Result collapse = SpectrumLevelCollapse.bestPerSpectrum(collapseInput);
		cb.onMessage("Spectrum-level collapse: " + collapseInput.size() + " candidates -> " + collapse.spectraCount
				+ " spectra (winner " + collapse.describeRanks() + ")");
		cb.onMessage(String.format("Classifier overruled the engine's top hit for %.1f%% of spectra",
				collapse.overruleFraction() * 100));

		// q-values recomputed on the collapsed set; the classifier's own figure
		// was per-candidate and is not comparable to a published PSM count.
		List<OpenSearchPsm> best = collapse.best;
		double[] bestScores = new double[best.size()];
		boolean[] bestTargets = new boolean[best.size()];
		for (int i = 0; i < best.size(); i++) {
			bestScores[i] = best.get(i).getComparableScore();
			bestTargets[i] = best.get(i).isTarget();
		}
		double[] bestQ = TargetDecoyFdr.qValues(bestScores, bestTargets);
		Map<OpenSearchPsm, Double> spectrumQValues = new IdentityHashMap<>();
		for (int i = 0; i < best.size(); i++) {
			best.get(i).setQValue(bestQ[i]);
			spectrumQValues.put(best.get(i), bestQ[i]);
		}
		int acceptedSpectra = TargetDecoyFdr.countAcceptedTargets(bestScores, bestTargets, 0.01);
		cb.onMessage("Classifier accepted " + acceptedSpectra + " target spectra at 1% FDR (per spectrum, "
				+ "not per candidate)");
		LOGGER.info("Per-spectrum FDR: " + acceptedSpectra + " target spectra accepted at 1% of "
				+ collapse.spectraCount + " spectra (scoring: " + engineDiagnostic + ")");

		if (joint) {
			// The comparison arms, each captured before the next step re-walks the same PSM objects:
			// closed only, open only, the cascade. The joint identification then scores everything
			// again on its own model and is what the workflow reports.
			String engines = params.getSearchEngines();
			IdentificationArms arms = new IdentificationArms(OpenSearchLayout.identificationArms(resultRoot, engines),
					OpenSearchLayout.identificationArmPeptides(resultRoot, engines),
					OpenSearchLayout.identificationArmsSummary(resultRoot, engines));
			arms.captureClosedOnly(closed, CASCADE_FDR);
			ModificationValidator.Result openOnly = ModificationValidator.validate(best, params.toValidatorConfig(),
					gate);
			arms.captureOpenOnly(best, openOnly, gate);
			ModificationValidator.Result cascadeResult = validateOpenStage(params, psms, best, closed, gate,
					substitutions, resultRoot, cb);
			arms.captureCascade(closed, best, cascadeResult, gate);
			ModificationValidator.Result validation = runJoint(params, psms, closedPsms, gate, spectraFiles,
					resultRoot, classifierConfig, cb, arms);
			alignSubstitutionVerdicts(params, validation, peptideIndex, cb);
			writePotentialModifications(fit, substitutions, validation,
					OpenSearchLayout.potentialModifications(resultRoot, params.getSearchEngines()));
			return validation;
		}
		ModificationValidator.Result validation = cascade
				? validateOpenStage(params, psms, best, closed, gate, substitutions, resultRoot, cb)
				: ModificationValidator.validate(best, params.toValidatorConfig(), gate);
		writeScoredPsms(psms, best, spectrumQValues, validation, closed,
				OpenSearchLayout.scoredPsms(resultRoot, params.getSearchEngines()));
		writePotentialModifications(fit, substitutions, validation,
				OpenSearchLayout.potentialModifications(resultRoot, params.getSearchEngines()));
		return validation;
	}

	/**
	 * DIA-NN predicted-spectrum and retention-time features for one search's PSMs, every reported
	 * rank, as external features for the cascade's classifiers.
	 *
	 * @return PSMs given an entropy; 0 when DIA-NN or the spectra are unavailable (warned)
	 */
	private static int attachCascadePredictions(DDAMagOpenSearchParameters params, List<OpenSearchPsm> psms,
			List<File> spectraFiles, File resultRoot, String search, ProgressCallback cb) {
		PredictedSpectrumFeatures.Inputs ms2 = new PredictedSpectrumFeatures.Inputs();
		ms2.diann = params.getDiannExecutable();
		ms2.mzmlBySample = VariantFormRanking.mzmlBySample(spectraFiles);
		ms2.threads = params.getThreadCount();
		ms2.fragmentTolerancePpm = params.getJointFragmentTolerancePpm();
		ms2.maxOpenRank = Integer.MAX_VALUE;
		ms2.workDir = new File(OpenSearchLayout.jointPredictionDir(resultRoot, params.getSearchEngines()),
				"cascade_" + search);
		String label = "Cascade " + search + " predictions";
		int n = PredictedSpectrumFeatures.attach(psms, ms2, label, cb::onMessage,
				message -> cb.onMessage("WARNING: " + message));
		cb.onMessage(String.format(java.util.Locale.ROOT, "%s: predicted-spectrum features for %,d of %,d PSMs", label,
				n, psms.size()));
		return n;
	}

	/** Spectrum- and peptide-level FDR of the closed-then-open identification. */
	private static final double CASCADE_FDR = 0.01;

	/**
	 * The joint closed+open identification ({@link JointIdentification}): writes its summary, its
	 * feature report and its scored table, and the comparison arms when {@code arms} is given.
	 *
	 * @param psms the open candidates with peak context and the peak filter applied
	 */
	/** Score features the v2 tree model may only reward (a better match never scores lower). */
	static final java.util.Set<String> V2_MONOTONE_FEATURES = java.util.Set.of("closed_hyperscore",
			"closed_neg_log10_expect", "open_hyperscore", "open_neg_log10_expect", "sage_hyperscore",
			"sage_neg_log10_expect", "hyperscore_margin", "ms2_entropy", "ms2_matched_fraction");

	/**
	 * The identification v2 joint model (design C3): monotone gradient-boosted trees, the audit
	 * harness's settings, all cores.
	 */
	static SemiSupervisedPsmClassifier.ModelFactory identificationV2Model() {
		return new SemiSupervisedPsmClassifier.ModelFactory() {
			@Override
			public weka.classifiers.Classifier create() {
				XGBoostTreeClassifier model = new XGBoostTreeClassifier(V2_MONOTONE_FEATURES);
				model.threads = Runtime.getRuntime().availableProcessors();
				return model;
			}

			@Override
			public String describe() {
				return "XGBoost trees (depth 4, 200 rounds, eta 0.1, lambda 1, monotone score features)";
			}
		};
	}

	/**
	 * Joint identification only: every reported PSM's substitution verdict becomes the identification's
	 * own call ({@link SubstitutionReassigner#alignWithIdentification}), so the variant outputs and the
	 * side quantification read the decision that kept the PSM out of the main result. Off with the
	 * substitution analysis.
	 */
	private static void alignSubstitutionVerdicts(DDAMagOpenSearchParameters params,
			ModificationValidator.Result validation, TrypticPeptideIndex peptideIndex, ProgressCallback cb) {
		if (!params.isSubstitutionAnalysisEnabled()) {
			return;
		}
		SubstitutionReassigner.Alignment a = SubstitutionReassigner.alignWithIdentification(reportedPsms(validation),
				peptideIndex);
		cb.onMessage(String.format(java.util.Locale.ROOT,
				"Substitution verdicts follow the identification: %,d set-aside PSMs whose only explanation is a "
						+ "substitution (%,d with the substituted peptide in the searched database); %,d peak-level "
						+ "verdicts removed from PSMs the identification did not call a substitution",
				a.called, a.inDatabase, a.cleared));
	}

	private static ModificationValidator.Result runJoint(DDAMagOpenSearchParameters params, List<OpenSearchPsm> psms,
			List<OpenSearchPsm> closedPsms, ModificationValidator.Gate gate, List<File> spectraFiles, File resultRoot,
			SemiSupervisedPsmClassifier.Config classifierConfig, ProgressCallback cb, IdentificationArms arms)
			throws Exception {
		String engines = params.getSearchEngines();
		JointIdentification.Config cfg = new JointIdentification.Config();
		cfg.fdr = CASCADE_FDR;
		cfg.acceptUnexplained = params.isJointAcceptUnexplained();
		cfg.withheldFeatures = params.jointWithheldFeatureSet();
		cfg.classifier.folds = classifierConfig.folds;
		cfg.classifier.maxIterations = classifierConfig.maxIterations;
		cfg.classifier.trainQValue = classifierConfig.trainQValue;
		if (params.isIdentificationV2()) {
			// docs/developer/identification-v2-design.md: precursor queries (C1), RESET (C2), trees (C3).
			cfg.precursorQueries = true;
			cfg.reset = true;
			cfg.classifier.modelFactory = identificationV2Model();
			cfg.classifier.unclampedLogits = true;
			// Section 7 of docs/developer/shift-localisation-explanation-design.md: localised-explained
			// unexplained-shift winners form their own FDR group; the model is unchanged.
			cfg.shiftLocaliser = new ShiftLocaliser(params.getJointFragmentTolerancePpm());
			cfg.localisationSpectra = VariantFormRanking.mzmlBySample(spectraFiles);
			cb.onMessage("Joint identification v2: precursor queries, RESET (decoy weight 2), "
					+ cfg.classifier.modelFactory.describe() + ", shift localisation (" + cfg.shiftLocaliser + ")");
		}
		PredictedSpectrumFeatures.Inputs ms2 = null;
		if (params.isJointMs2Features()) {
			ms2 = new PredictedSpectrumFeatures.Inputs();
			ms2.diann = params.getDiannExecutable();
			ms2.mzmlBySample = VariantFormRanking.mzmlBySample(spectraFiles);
			ms2.threads = params.getThreadCount();
			ms2.fragmentTolerancePpm = params.getJointFragmentTolerancePpm();
			ms2.maxOpenRank = params.getJointMs2MaxOpenRank();
			ms2.workDir = OpenSearchLayout.jointPredictionDir(resultRoot, engines);
		}
		cb.onMessage("Joint identification: closed and open candidates in one pool, one model, FDR per evidence group"
				+ (ms2 == null ? " (predicted-spectrum features off)" : ""));
		JointIdentification.Result joint = JointIdentification.run(psms, closedPsms, gate, cfg, ms2, cb::onMessage,
				message -> cb.onMessage("WARNING: " + message));

		File summary = OpenSearchLayout.jointSummary(resultRoot, engines);
		summary.getParentFile().mkdirs();
		try (PrintWriter writer = new PrintWriter(summary, "UTF-8")) {
			writer.println("# Joint closed+open identification (docs/developer/joint-closed-open-identification-plan.md).");
			writer.println("# Counts in [targets, decoys]. Peptides: base peptide with I = L, each in the highest-priority");
			writer.println("# group where it has an accepted PSM (exact > open_unmodified > open_accepted_shift).");
			writer.println("name\tvalue");
			for (String line : joint.summary) {
				writer.println(line.contains("\t") ? line : "model\t" + line);
			}
		}
		JointIdentification.writeFeatureReport(joint, OpenSearchLayout.jointFeatureReport(resultRoot, engines));
		if (!joint.peakFragmentTests.isEmpty()) {
			PeakFragmentTest.write(joint.peakFragmentTests, OpenSearchLayout.peakFragmentTest(resultRoot, engines));
		}
		writeJointScoredPsms(joint, OpenSearchLayout.scoredPsms(resultRoot, engines));
		if (arms != null) {
			arms.captureJoint(joint);
			arms.finish();
			cb.onMessage("Joint identification: comparison arms written to "
					+ OpenSearchLayout.identificationArmsSummary(resultRoot, engines).getName());
		}
		cb.onMessage(String.format(java.util.Locale.ROOT,
				"Joint identification: %,d identifications (%,d exact, %,d open) at %.0f%% FDR per group",
				joint.validation.validated.size(),
				joint.validation.validated.size() - joint.validation.peakValidated.size(),
				joint.validation.peakValidated.size(), CASCADE_FDR * 100));
		return joint.validation;
	}

	/**
	 * Audit table of the joint identification: one row per candidate (spectrum, peptide), with the
	 * searches that proposed it, each search's own scores, the joint score, the winner flag, the
	 * group and both q-values. {@link OpenSearchGenomeEvidence} reads its winners.
	 */
	private static void writeJointScoredPsms(JointIdentification.Result joint, File out) throws IOException {
		out.getParentFile().mkdirs();
		String[] featureNames = joint.features.featureNames();
		List<String> ms2Columns = new ArrayList<>();
		for (String name : featureNames) {
			if (name.startsWith("ms2_") || name.startsWith("rt_")) {
				ms2Columns.add(name);
			}
		}
		int[] ms2Index = new int[ms2Columns.size()];
		for (int k = 0; k < ms2Index.length; k++) {
			ms2Index[k] = Arrays.asList(featureNames).indexOf(ms2Columns.get(k));
		}
		try (PrintWriter writer = new PrintWriter(new java.io.BufferedWriter(
				new java.io.OutputStreamWriter(new java.io.FileOutputStream(out), java.nio.charset.StandardCharsets.UTF_8),
				1 << 20))) {
			StringBuilder header = new StringBuilder("file\tscan\tpeptide\tprotein\tis_target\tcharge\t"
					+ "precursor_neutral_mass\tretention_time_sec\tfound_by\tid_group\tclosed_hit_rank\t"
					+ "open_hit_rank\tsage_hit_rank\tmass_diff\tnormalized_mass_diff\tpeak_id\t"
					+ "peak_hypothesis_accepted\tclosed_hyperscore\tclosed_expect\topen_hyperscore\topen_expect\t"
					+ "delta_mass_score_gain\tlocalized_sequence\tclass_score\tclass_logit\tis_spectrum_winner\t"
					+ "spectrum_q_value\tpeptide_q_value\tvalidated\tstage");
			boolean v2Columns = joint.trainingDecoy != null;
			for (String name : ms2Columns) {
				header.append('\t').append(name);
			}
			if (v2Columns) {
				header.append("\tquery_reassigned\treset_role\tshift_state\tshift_site\tshift_site_probability"
						+ "\tshift_support\tshift_explanation\tshift_class");
			}
			writer.println(header);
			for (int i = 0; i < joint.candidates.size(); i++) {
				JointCandidate c = joint.candidates.get(i);
				OpenSearchPsm rep = c.representative();
				OpenSearchPsm open = c.openPrimary();
				StringBuilder row = new StringBuilder(256);
				row.append(c.getFileName()).append('\t').append(c.getScan()).append('\t').append(c.getPeptide())
						.append('\t').append(nullSafe(rep.getProtein())).append('\t').append(c.isTarget())
						.append('\t').append(rep.getCharge()).append('\t').append(rep.getPrecursorNeutralMass())
						.append('\t').append(rep.getRetentionTime()).append('\t').append(c.support()).append('\t')
						.append(c.getGroup()).append('\t').append(c.closed() == null ? "" : c.closed().getHitRank())
						.append('\t').append(c.fragger() == null ? "" : c.fragger().getHitRank()).append('\t')
						.append(c.sage() == null ? "" : c.sage().getHitRank()).append('\t').append(rep.getMassDiff())
						.append('\t').append(rep.getNormalizedMassDiff()).append('\t')
						.append(open == null ? -1 : open.getPeakId()).append('\t')
						.append(open != null && open.isInValidatedPeak()).append('\t')
						.append(c.closed() == null ? "" : c.closed().getHyperscore()).append('\t')
						.append(c.closed() == null ? "" : c.closed().getExpect()).append('\t')
						.append(c.fragger() == null ? "" : c.fragger().getHyperscore()).append('\t')
						.append(c.fragger() == null ? "" : c.fragger().getExpect()).append('\t')
						.append(open == null ? "" : open.getDeltaMassScoreGain()).append('\t')
						.append(open == null ? "" : nullSafe(open.getLocalizedSequence())).append('\t')
						.append(joint.scores[i]).append('\t').append(joint.logits[i]).append('\t')
						.append(joint.winner[i]).append('\t').append(valueOrBlank(nanToNull(joint.psmQ[i])))
						.append('\t').append(valueOrBlank(nanToNull(joint.peptideQ[i]))).append('\t')
						.append(joint.kept[i]).append('\t').append(rep.isClosedStage() ? "closed" : "open");
				if (ms2Index.length > 0) {
					double[] v = joint.features.values(i);
					for (int k : ms2Index) {
						row.append('\t').append(Double.isNaN(v[k]) ? "" : String.valueOf(v[k]));
					}
				}
				if (v2Columns) {
					row.append('\t').append(c.isReassigned()).append('\t').append(c.isTarget() ? "target"
							: joint.trainingDecoy[i] ? "training_decoy" : "estimating_decoy");
					// Shift localisation of unexplained-shift winners; the site is 1-based, as residues are named.
					ShiftLocaliser.Result loc = c.getLocalisation();
					if (loc == null) {
						row.append("\t\t\t\t\t\t");
					} else {
						row.append('\t').append(loc.state.name().toLowerCase(java.util.Locale.ROOT)).append('\t')
								.append(loc.site < 0 ? "" : String.valueOf(loc.site + 1)).append('\t')
								.append(valueOrBlank(nanToNull(loc.siteProbability))).append('\t')
								.append(valueOrBlank(nanToNull(loc.support))).append('\t')
								.append(loc.explanationText()).append('\t').append(loc.explanationClassText());
					}
				}
				writer.println(row);
			}
		}
	}

	private static Double nanToNull(double v) {
		return Double.isNaN(v) ? null : v;
	}

	/**
	 * Open stage and peptide filter of the closed-then-open identification: the open winners of
	 * spectra the closed stage did not accept go through {@link ModificationValidator#validateRemaining};
	 * the closed and open identifications then go through the per-class peptide walks and the
	 * union check ({@link ClosedOpenCascade#filterPeptides}). Writes the run summary.
	 */
	private static ModificationValidator.Result validateOpenStage(DDAMagOpenSearchParameters params,
			List<OpenSearchPsm> candidates, List<OpenSearchPsm> best, ClosedOpenCascade.ClosedStage closed,
			ModificationValidator.Gate gate, SubstitutionReassigner.Result substitutions,
			File resultRoot, ProgressCallback cb) throws IOException {
		List<OpenSearchPsm> remaining = new ArrayList<>();
		int closedTargetSpectra = 0;
		int openWinnerDiffers = 0;
		for (OpenSearchPsm psm : best) {
			OpenSearchPsm closedPsm = closed.acceptedBySpectrum.get(ClosedOpenCascade.spectrumKey(psm));
			if (closedPsm == null) {
				remaining.add(psm);
				continue;
			}
			if (closedPsm.isTarget()) {
				closedTargetSpectra++;
				if (!psm.isTarget() || !ClosedOpenCascade.basePeptide(psm.getPeptide())
						.equals(ClosedOpenCascade.basePeptide(closedPsm.getPeptide()))) {
					openWinnerDiffers++;
				}
			}
		}
		ModificationValidator.Result open = ModificationValidator.validateRemaining(remaining,
				params.toValidatorConfig(), gate);

		Map<Integer, ModificationValidator.MassClass> classes = ModificationValidator.massClasses(gate);
		Map<ClosedOpenCascade.IdClass, List<OpenSearchPsm>> pools = new EnumMap<>(ClosedOpenCascade.IdClass.class);
		Map<ClosedOpenCascade.IdClass, List<OpenSearchPsm>> accepted = new EnumMap<>(
				ClosedOpenCascade.IdClass.class);
		for (ClosedOpenCascade.IdClass c : ClosedOpenCascade.IdClass.values()) {
			pools.put(c, new ArrayList<>());
			accepted.put(c, new ArrayList<>());
		}
		pools.get(ClosedOpenCascade.IdClass.CLOSED).addAll(closed.winners);
		accepted.get(ClosedOpenCascade.IdClass.CLOSED).addAll(closed.accepted);
		Set<OpenSearchPsm> openAccepted = Collections.newSetFromMap(new IdentityHashMap<>());
		openAccepted.addAll(open.peakValidated);
		for (OpenSearchPsm psm : remaining) {
			if (!gate.acceptedPeakIds.contains(psm.getPeakId())) {
				continue;
			}
			ModificationValidator.MassClass massClass = classes.get(psm.getPeakId());
			ClosedOpenCascade.IdClass c = massClass == ModificationValidator.MassClass.ZERO
					? ClosedOpenCascade.IdClass.OPEN_UNMODIFIED
					: massClass == ModificationValidator.MassClass.ESTABLISHED ? ClosedOpenCascade.IdClass.OPEN_ACCEPTED_SHIFT
							: null;
			if (c == null) {
				continue;
			}
			pools.get(c).add(psm);
			if (openAccepted.contains(psm)) {
				accepted.get(c).add(psm);
			}
		}
		ClosedOpenCascade.PeptideFilter peptides = ClosedOpenCascade.filterPeptides(pools, accepted, CASCADE_FDR);
		List<OpenSearchPsm> validated = new ArrayList<>(peptides.kept);
		List<OpenSearchPsm> openKept = new ArrayList<>();
		int closedKept = 0;
		for (OpenSearchPsm psm : validated) {
			if (psm.isClosedStage()) {
				closedKept++;
			} else {
				openKept.add(psm);
			}
		}
		cb.onMessage(String.format(java.util.Locale.ROOT,
				"Closed then open: %,d closed + %,d open identifications after the peptide filter "
						+ "(peptide threshold factor %.3f, peptide union FDR %.2f%%); %,d unexplained-shift "
						+ "PSMs reported, not identified; the open winner differs on %,d of %,d closed target spectra",
				closedKept, openKept.size(), peptides.factor, 100 * peptides.unionFdrAfter(),
				open.unexplainedIdentifications.size(), openWinnerDiffers, closedTargetSpectra));

		File summary = OpenSearchLayout.closedOpenSummary(resultRoot, params.getSearchEngines());
		summary.getParentFile().mkdirs();
		try (PrintWriter writer = new PrintWriter(summary, "UTF-8")) {
			writer.println("# Closed-then-open identification (docs/developer/closed-then-open-identification-plan.md).");
			writer.println("# Counts in [targets, decoys] are per class. Peptides: base peptide with I = L.");
			writer.println("name\tvalue");
			writer.println("closed_spectra\t" + closed.winners.size());
			writer.println("closed_spectra_accepted\t[" + closed.acceptedTargets + ", " + closed.acceptedDecoys + "]");
			writer.println("open_spectra\t" + best.size());
			writer.println("open_spectra_remaining\t" + remaining.size());
			writer.println("closed_target_spectra_with_differing_open_winner\t" + openWinnerDiffers + "/"
					+ closedTargetSpectra);
			writer.println("open_unexplained_reported\t" + open.unexplainedIdentifications.size());
			for (String line : ClosedOpenCascade.describe(peptides)) {
				writer.println(line);
			}
			writer.println("validated_closed\t" + closedKept);
			writer.println("validated_open\t" + openKept.size());
			writeSubstitutionAgreement(writer, candidates, closed, substitutions);
		}
		return new ModificationValidator.Result(validated, openKept, open.unexplainedIdentifications, open.verdicts,
				open.diagnostic);
	}

	/**
	 * Per substitution-shaped peak: its rank-1 targets on a spectrum the closed stage accepted as a
	 * target, agreeing (same base peptide) or contradicted. Reported only; decides nothing.
	 */
	private static void writeSubstitutionAgreement(PrintWriter writer, List<OpenSearchPsm> candidates,
			ClosedOpenCascade.ClosedStage closed, SubstitutionReassigner.Result substitutions) {
		Map<Integer, int[]> counts = new TreeMap<>();
		for (SubstitutionReassigner.PeakSummary peak : substitutions.peaks) {
			if (peak.substitutionShaped) {
				counts.put(peak.peakId, new int[2]);
			}
		}
		for (OpenSearchPsm psm : candidates) {
			int[] c = counts.get(psm.getPeakId());
			if (c == null || psm.getHitRank() != 1 || !psm.isTarget()) {
				continue;
			}
			OpenSearchPsm closedPsm = closed.acceptedBySpectrum.get(ClosedOpenCascade.spectrumKey(psm));
			if (closedPsm == null || !closedPsm.isTarget()) {
				continue;
			}
			boolean same = ClosedOpenCascade.basePeptide(psm.getPeptide())
					.equals(ClosedOpenCascade.basePeptide(closedPsm.getPeptide()));
			c[same ? 0 : 1]++;
		}
		for (SubstitutionReassigner.PeakSummary peak : substitutions.peaks) {
			int[] c = counts.get(peak.peakId);
			if (c == null) {
				continue;
			}
			writer.println(String.format(java.util.Locale.ROOT,
					"substitution_peak_%d\tcentroid=%+.4f closed_agree=%d closed_contradicted=%d",
					peak.peakId, peak.centroid, c[0], c[1]));
		}
	}

	/**
	 * Refines one stage's pepXMLs with Crystal-C, and returns the files the rest
	 * of the workflow should read.
	 *
	 * <h2>Why</h2>
	 * An open search reports {@code massdiff = precursor_neutral_mass -
	 * calc_neutral_pep_mass}. When the instrument picked the wrong peak of the
	 * precursor's isotope cluster, or a co-eluting neighbour inside the isolation
	 * window, that difference is a precursor-assignment error a few Da from zero
	 * - not a modification. Crystal-C re-reads the MS1 scan, re-picks the
	 * monoisotopic peak and rewrites {@code precursor_neutral_mass}.
	 *
	 * <h2>What it is worth, measured</h2>
	 * On 12mix (2026-09-14), over the nine stage-2 pepXMLs and their Crystal-C
	 * outputs, rank-1 PSMs at 1% FDR on hyperscore:
	 * <pre>
	 *   population            before              after
	 *   |delta| &lt; 0.05 Da     310,049 / 96.98%    581,451 / 97.03%   (+87.5%)
	 *   |delta| &lt;= 5 Da       659,561 / 96.98%    675,149 / 96.98%   (+2.4%)
	 *   all rank-1            577,440 / 96.35%    577,440 / 96.35%   (identical)
	 * </pre>
	 * (accepted target PSMs / share whose protein set includes one of the 12
	 * spiked genomes.) Read the last line first: <b>Crystal-C creates no
	 * identifications.</b> It changes no score and, at rank 1, essentially no
	 * peptide. What it does is move 271,402 correct PSMs and 16,644 peptides out
	 * of spurious delta-mass peaks and into the zero peak, at unchanged PSM
	 * purity. Everything this step buys comes from that relocation - which is why
	 * it matters most at stage 1, whose admission rule is the zero family plus
	 * tested UNIMOD shifts, and nothing else.
	 * <p>
	 * It is not free: decoys inside the rank-1 zero window rose from 1,192 to
	 * 6,465 (0.38% to 1.10% of the window), so roughly 1.9% of what it relocates
	 * is decoy and the window's internal FDR threshold tightens accordingly.
	 *
	 * <h2>Ranks above 1</h2>
	 * Crystal-C rewrites the {@code massdiff} attribute on <b>rank 1 only</b>.
	 * Measured on {@code YS-194-1-1_c.pepXML}: rank 1 agrees with
	 * {@code precursor_neutral_mass - calc_neutral_pep_mass} to 1e-5 Da, while
	 * 29,662 of 137,739 rank-2 hits disagree by more than 0.005 Da, up to 499 Da,
	 * and ranks 3-5 likewise. The attribute is stale on the lower ranks because
	 * the spectrum's precursor moved underneath it. This is why
	 * {@code OpenSearchPepxmlReader} derives the delta rather than reading the
	 * attribute - without that, stage 2 would score corrected precursors against
	 * uncorrected deltas for about a fifth of its candidates.
	 *
	 * <h2>Ordering</h2>
	 * Stage 1: immediately after the marker search and before MSBooster, so
	 * genome selection sees corrected masses. MSBooster is unaffected - it reads
	 * the {@code .pin}, which Crystal-C does not touch, and the rescorer joins it
	 * back by scan.
	 * <p>
	 * Stage 2: before Sage and before anything reads a delta mass. Sage is
	 * unaffected: Crystal-C takes pepXML only, so a MSFragger+Sage run refines
	 * the MSFragger half and leaves Sage's precursors as searched - a known
	 * asymmetry, not a silent one.
	 *
	 * <h2>Resume</h2>
	 * The corrected pepXML is this step's sentinel, gated on the search output it
	 * was derived from, so a re-run of the search invalidates it. When the step is
	 * switched off the uncorrected pepXMLs are returned unchanged and any stale
	 * {@code _c.pepXML} beside them is ignored rather than picked up.
	 *
	 * @param searchDir  {@link OpenSearchLayout#MARKER_SEARCH_DIR} or
	 *                   {@link OpenSearchLayout#SAMPLE_DB_SEARCH_DIR}; picks both
	 *                   the folder and which database was searched.
	 * @param markerTd   the stage-1 target-decoy marker database, or null at
	 *                   stage 2, where each sample's refined FASTA is used.
	 */
	private static List<File> runCrystalC(DDAMagOpenSearchParameters params, File resultRoot,
			List<String[]> experiments, List<File> pepXmls, String searchDir, File markerTd, ProgressCallback cb)
			throws IOException {
		boolean marker = OpenSearchLayout.MARKER_SEARCH_DIR.equals(searchDir);
		if (!params.isCrystalcEnabled()) {
			cb.onMessage("Crystal-C precursor refinement is switched off for " + searchDir + "; delta masses will "
					+ "include isotope-selection and co-isolation artifacts.");
			return pepXmls;
		}
		File jar = params.getCrystalcJar();
		if (jar == null || !jar.isFile()) {
			throw new IOException("Crystal-C precursor refinement is enabled but its jar was not found. Set "
					+ "Resources > FragPipe to an installation whose tools folder contains "
					+ "original-crystalc-*.jar (with batmass-io and grppr beside it), or switch the step off "
					+ "deliberately - running without it reports precursor-assignment artifacts as "
					+ "modifications.");
		}
		File spectraDir = readableSpectraDir(experiments, resultRoot);
		if (spectraDir == null || !spectraDir.isDirectory()) {
			throw new IOException("Crystal-C re-reads the MS1 scans, so it needs the folder holding the spectra; "
					+ "resolved to " + spectraDir + ".");
		}

		List<File> refined = new ArrayList<>(experiments.size());
		for (int i = 0; i < experiments.size(); i++) {
			String sample = baseName(experiments.get(i)[0]);
			File pepXml = pepXmls.get(i);
			File workDir = marker ? OpenSearchLayout.markerDir(resultRoot, sample)
					: OpenSearchLayout.sampleDbDir(resultRoot, sample);
			File corrected = marker ? OpenSearchLayout.markerCorrectedPepXml(resultRoot, sample)
					: OpenSearchLayout.sampleDbCorrectedPepXml(resultRoot, sample);
			if (OpenSearchLayout.isSatisfied(corrected, pepXml)) {
				cb.onMessage("Crystal-C output already exists, skipping: " + sample + "/" + searchDir);
				refined.add(corrected);
				continue;
			}
			if (OpenSearchLayout.isSatisfied(corrected)) {
				cb.onMessage(sample + "/" + searchDir
						+ ": Crystal-C output is older than its search; re-running it");
			}

			CrystalCHeadlessRunner.Settings settings = new CrystalCHeadlessRunner.Settings();
			settings.rawFileLocation = spectraDir;
			settings.rawFileExtension = extensionOf(
					OpenSearchSpectra.readable(new File(experiments.get(i)[0]), resultRoot).getName());
			settings.outputLocation = workDir;
			// The database this search ran against. Crystal-C needs the sequences
			// for its missed-cleavage artifact class and refuses an empty value
			// outright. At stage 1 that is the shared marker target-decoy FASTA,
			// which is MAG-scale - it is loaded whole, so this is the step's
			// memory floor, not the spectra.
			settings.fasta = marker ? markerTd : OpenSearchLayout.refinedFasta(resultRoot, sample);
			settings.threads = Math.max(1, params.getThreadCount());
			settings.minCharge = params.getCrystalcMinCharge();
			settings.maxCharge = params.getCrystalcMaxCharge();
			settings.isotopeNumber = params.getCrystalcIsotopeNumber();
			settings.precursorMassPpm = params.getCrystalcPrecursorMassPpm();
			settings.precursorIsolationWindow = params.getCrystalcIsolationWindow();

			File crystalcParams = CrystalCHeadlessRunner.writeParams(workDir, pepXml, settings);
			cb.onMessage("Crystal-C: refining precursor masses for " + sample + "/" + searchDir);
			int exit;
			try {
				exit = CrystalCHeadlessRunner.run(jar, params.getCrystalcMaxHeap(), crystalcParams, pepXml, cb);
			} catch (InterruptedException e) {
				Thread.currentThread().interrupt();
				throw new IOException("Interrupted while running Crystal-C", e);
			}
			if (exit != 0) {
				throw new IOException("Crystal-C exited " + exit + " for " + sample + "/" + searchDir
						+ ". Its parameters are in " + crystalcParams.getAbsolutePath() + "; it calls "
						+ "System.exit(1) from its own parameter parser, so check those first.");
			}
			if (!OpenSearchLayout.isSatisfied(corrected)) {
				throw new IOException("Crystal-C reported success but wrote no " + corrected.getName() + " for "
						+ sample + ". It skips any input already ending _c.pepXML, and it reads the spectra from "
						+ spectraDir.getAbsolutePath() + " - check its output above.");
			}
			refined.add(corrected);
		}
		cb.onMessage("Crystal-C: " + refined.size() + " refined " + searchDir + " pepXML(s); everything "
				+ "downstream reads these, not the raw search output.");
		return refined;
	}

	/**
	 * Peptide space of the databases this run actually searched, for deciding
	 * whether a proposed substituted sequence was ever a candidate.
	 * <p>
	 * Built from the per-sample refined FASTAs - the same files stage 2 searched -
	 * so a "not in the database" verdict means exactly that and not "not in some
	 * other database". Returns null when the analysis is switched off or no
	 * refined database can be found, which {@code SubstitutionReassigner} reads as
	 * "cannot test membership" and reports every candidate as absent rather than
	 * silently claiming it is present.
	 */
	private static TrypticPeptideIndex buildPeptideIndex(DDAMagOpenSearchParameters params, File resultRoot,
			ProgressCallback cb) {
		if (!params.isSubstitutionAnalysisEnabled()) {
			return null;
		}
		List<File> fastas = refinedFastas(params, resultRoot);
		if (fastas.isEmpty()) {
			cb.onMessage("Substitution analysis: no refined database found, so candidate sequences cannot be "
					+ "tested for membership; all candidates will be reported as absent from the search space.");
			return null;
		}
		TrypticPeptideIndex.Config cfg = new TrypticPeptideIndex.Config();
		cfg.missedCleavages = params.getMissCleavages();
		cfg.decoyTag = TargetDecoyFasta.DECOY_TAG;
		try {
			long start = System.currentTimeMillis();
			TrypticPeptideIndex index = TrypticPeptideIndex.build(fastas, cfg);
			cb.onMessage(String.format(java.util.Locale.ROOT,
					"Substitution analysis: indexed %,d distinct tryptic peptides from %,d proteins across %d "
							+ "refined database(s) in %.1f s",
					index.size(), index.proteinCount(), fastas.size(),
					(System.currentTimeMillis() - start) / 1000.0));
			return index;
		} catch (IOException e) {
			LOGGER.warn("Could not index the refined databases for substitution analysis", e);
			cb.onMessage("Substitution analysis: could not index the refined databases (" + e.getMessage()
					+ "); candidates will be reported as absent from the search space.");
			return null;
		}
	}

	/** The per-sample refined FASTAs stage 2 searched, those that exist on disk. */
	private static List<File> refinedFastas(DDAMagOpenSearchParameters params, File resultRoot) {
		List<File> fastas = new ArrayList<>();
		if (params.getExperimentData() == null) {
			return fastas;
		}
		for (String[] experiment : params.getExperimentData()) {
			File refined = OpenSearchLayout.refinedFasta(resultRoot, baseName(experiment[0]));
			if (refined.isFile()) {
				fastas.add(refined);
			}
		}
		return fastas;
	}

	/**
	 * Calls sequence variants and writes {@code sequence_variants.tsv} and
	 * {@code sequence_variants.fasta}. Stale copies from an earlier run are
	 * removed when the call is switched off or yields nothing, so a file on disk
	 * always belongs to the run beside it. A failure here is logged and reported,
	 * never fatal: the validated PSMs do not depend on it.
	 */
	/**
	 * The identified PSMs exported as the unmodified peptide, closed and open alike: what the variant
	 * caller reads to say a reference peptide is also observed. An accepted shift exported unmodified
	 * for want of fragment support is left out, since its precursor does carry a shift. Null unless
	 * the result separates identifications from the variant evidence (the joint identification), so
	 * every other identification keeps the caller's own reading.
	 */
	private static List<OpenSearchPsm> unmodifiedIdentifications(ModificationValidator.Result validation) {
		if (validation.diagnostic == validation.validated) {
			return null;
		}
		List<OpenSearchPsm> out = new ArrayList<>();
		for (OpenSearchPsm psm : validation.validated) {
			if (psm != null && psm.getPeptide() != null && !psm.isShiftFoldedUnmodified()
					&& psm.getPeptide().equals(OpenSearchPostAnalysisTsvExporter.canonicalModifiedSequenceOf(psm))) {
				out.add(psm);
			}
		}
		return out;
	}

	private static SequenceVariantCaller.Result callSequenceVariants(DDAMagOpenSearchParameters params,
			File resultRoot, List<OpenSearchPsm> validated, List<OpenSearchPsm> references, ProgressCallback cb) {
		File table = OpenSearchLayout.sequenceVariants(resultRoot, params.getSearchEngines());
		File fasta = OpenSearchLayout.sequenceVariantsFasta(resultRoot, params.getSearchEngines());
		SequenceVariantCaller.Config cfg = params.toSequenceVariantConfig();
		cfg.decoyTag = TargetDecoyFasta.DECOY_TAG;
		try {
			List<File> fastas = refinedFastas(params, resultRoot);
			SequenceVariantCaller.Result result = SequenceVariantCaller.call(validated, references, fastas, cfg,
					cb::onMessage);
			if (!cfg.enabled || result.variants.isEmpty()) {
				Files.deleteIfExists(table.toPath());
				Files.deleteIfExists(fasta.toPath());
				return result;
			}
			SequenceVariantCaller.writeVariantTable(result, cfg, table);
			int proteins = SequenceVariantCaller.writeVariantFasta(result, fastas, cfg, fasta);
			cb.onMessage("Sequence variants: " + result.variants.size() + " site(s) written to " + table.getName()
					+ "; " + proteins + " variant protein(s) for the recurrent classes written to " + fasta.getName());
			return result;
		} catch (IOException | RuntimeException e) {
			LOGGER.warn("Sequence-variant calling failed; validated PSMs are unaffected", e);
			cb.onMessage("WARNING: sequence-variant calling failed (" + e.getMessage()
					+ "); validated PSMs are unaffected.");
			return SequenceVariantCaller.Result.empty();
		}
	}

	/** Extension of a spectrum file, without the dot; Crystal-C needs it by name. */
	private static String extensionOf(String path) {
		String name = new File(path).getName();
		int dot = name.lastIndexOf('.');
		return dot < 0 ? "mzML" : name.substring(dot + 1);
	}

	/**
	 * Runs Sage against each sample's refined database, when the engine mode asks
	 * for it. Returns the per-sample {@code results.sage.tsv} files, or an empty
	 * list for a MSFragger-only run.
	 * <p>
	 * Sage searches the <b>same</b> refined FASTA MSFragger used, so the two
	 * engines differ only in how they score the same candidate space. Output goes
	 * to {@code sample_db_sage/}, which leaves MSFragger's {@code sample_db/}
	 * untouched and reusable — that is what lets a second run add Sage without
	 * repeating round one or MSFragger's round two.
	 */
	/** Whole-file read, used only to compare a previous Sage config with the new one. */
	private static String readFileText(File file) {
		try {
			return new String(Files.readAllBytes(file.toPath()), java.nio.charset.StandardCharsets.UTF_8);
		} catch (IOException e) {
			return "";
		}
	}

	private static List<File> runSageSampleDbSearches(DDAMagOpenSearchParameters params, File resultRoot,
			List<String[]> experiments, ProgressCallback cb) throws Exception {
		if (!params.usesSage()) {
			return Collections.emptyList();
		}
		File sageExe = params.getSageExecutable() == null ? null : new File(params.getSageExecutable());
		if (sageExe == null || !sageExe.exists()) {
			throw new IOException("Sage executable not found: "
					+ (sageExe == null ? "not configured" : sageExe.getAbsolutePath())
					+ ". Set it under Resources, or select MSFragger only.");
		}
		List<File> results = new ArrayList<>();
		for (String[] experiment : experiments) {
			String sample = baseName(experiment[0]);
			File dir = OpenSearchLayout.sampleDbDir(resultRoot, sample, "sage");
			File sentinel = SageOpenSearchRunner.resultsFile(dir);
			File refined = OpenSearchLayout.refinedFasta(resultRoot, sample);
			if (!refined.exists()) {
				throw new IOException("Refined database missing for " + sample + ": " + refined.getAbsolutePath());
			}
			dir.mkdirs();

			SageOpenSearchRunner.Settings settings = new SageOpenSearchRunner.Settings();
			settings.targetDecoyFasta = refined;
			settings.precursorMassLower = params.getPrecursorMassLower();
			settings.precursorMassUpper = params.getPrecursorMassUpper();
			settings.missedCleavages = params.getMissCleavages();
			// Mirrors MSFragger's round-2 topN (not the general outputReportTopN):
			// both are the same round-two per-sample search, just different engines,
			// so they should report the same number of candidates per spectrum.
			// Found 2026-09-18 during the round2OutputReportTopN UI-control fix:
			// this was still reading the general field, so it never tracked round1/
			// round2 changes at all.
			settings.reportPsms = params.getRound2OutputReportTopN();
			settings.decoyTag = TargetDecoyFasta.DECOY_TAG;
			settings.threads = params.getThreadCount();
			String[] enzymeInfo = MSFraggerHeadlessRunner.enzymeInfo(params.getEnzyme());
			settings.enzymeCleaveAt = enzymeInfo[1];
			settings.enzymeRestrict = enzymeInfo[2] == null ? "" : enzymeInfo[2];

			// The stage-2 MSFragger search for this sample has already run and left
			// its fragger.params beside its pepXML. Everything that defines the
			// search space - the open window, fixed modifications, peptide length
			// and mass limits, fragment tolerance, reported candidates - is copied
			// from that file rather than restated here, so an edit to the MSFragger
			// template cannot leave Sage searching a different problem. In
			// particular this is where the fixed carbamidomethyl comes from:
			// without it Sage searches an unmodified cysteine backbone.
			File fraggerParams = OpenSearchLayout.sampleDbParams(resultRoot, sample);
			if (!fraggerParams.exists()) {
				// Not a warning: the whole point of the second engine is that it
				// searches the same problem, and without this file that cannot be
				// guaranteed. Falling back to defaults would produce results that
				// look fine and are not comparable.
				throw new IOException("MSFragger parameters missing for " + sample + ": "
						+ fraggerParams.getAbsolutePath() + ". Sage mirrors its search space from this file, so the"
						+ " stage-2 MSFragger search must have run for this sample first.");
			}
			SageOpenSearchRunner.applyMsfraggerParams(settings, fraggerParams);

			// Staleness, by configuration rather than by timestamp. A changed open
			// window or fixed modification leaves the previous results.sage.tsv
			// intact and newer than everything that produced it, so only comparing
			// the config text catches it. This is what makes a settings fix
			// actually reach the results on a re-run instead of being skipped.
			File configFile = new File(dir, SageOpenSearchRunner.CONFIG_FILE);
			String wanted = SageOpenSearchRunner.buildConfigJson(dir, settings,
					Collections.singletonList(new File(experiment[0])));
			boolean configUnchanged = configFile.exists()
					&& wanted.trim().equals(readFileText(configFile).trim());
			if (configUnchanged && OpenSearchLayout.isSatisfied(sentinel, refined)) {
				cb.onMessage("Sage search already exists under the same configuration, skipping: " + sample);
				results.add(sentinel);
				continue;
			}
			if (OpenSearchLayout.isSatisfied(sentinel)) {
				cb.onMessage(sample + ": " + (configUnchanged
						? "Sage results are older than the refined database"
						: "the Sage configuration has changed since these results were produced")
						+ "; re-running the search rather than reusing them");
			}

			File config = SageOpenSearchRunner.writeConfig(dir, settings,
					Collections.singletonList(new File(experiment[0])));
			cb.onMessage(String.format(
					"%s: Sage mirroring %s - open window %.1f to %.1f Da (written as da [%.1f, %.1f]),"
							+ " fixed mods %s, %d candidates per spectrum",
					sample, fraggerParams.getName(), settings.precursorMassLower, settings.precursorMassUpper,
					-settings.precursorMassUpper, -settings.precursorMassLower, settings.staticMods,
					settings.reportPsms));
			cb.onMessage(sample + ": running Sage against the refined database ...");
			int exit = SageOpenSearchRunner.run(sageExe, config, settings.threads, cb);
			if (exit != 0) {
				throw new IOException("Sage search failed for " + sample + " (exit=" + exit + ")");
			}
			if (!sentinel.exists() || sentinel.length() == 0) {
				throw new IOException("Missing Sage results after search: " + sentinel.getAbsolutePath());
			}
			results.add(sentinel);
		}
		return results;
	}

	private static MSFraggerHeadlessRunner.OpenSearchSettings baseSettings(DDAMagOpenSearchParameters params,
			File fasta) {
		MSFraggerHeadlessRunner.OpenSearchSettings settings = new MSFraggerHeadlessRunner.OpenSearchSettings();
		settings.targetDecoyFasta = fasta;
		settings.threads = params.getThreadCount();
		// Name AND cleavage rules, from the shared mapping: setting only the name
		// would leave a non-tryptic enzyme digesting at KR.
		String[] enzymeInfo = MSFraggerHeadlessRunner.enzymeInfo(params.getEnzyme());
		settings.enzymeName = enzymeInfo[0];
		settings.cutAfter = enzymeInfo[1];
		settings.notCutBefore = enzymeInfo[2];
		settings.missedCleavages = params.getMissCleavages();
		settings.precursorMassLower = params.getPrecursorMassLower();
		settings.precursorMassUpper = params.getPrecursorMassUpper();
		settings.outputReportTopN = params.getOutputReportTopN();
		settings.decoyTag = TargetDecoyFasta.DECOY_TAG;
		return settings;
	}

	/** Copies the GUI/CLI genome-stop fields onto the selector request. */
	static void copyGenomeSelection(MagGenomeSelectionService.Request request, DDAMagOpenSearchParameters params) {
		request.genomeSelectorMode = params.getGenomeSelectorMode();
		request.genomeSelectionStrategy = params.getGenomeSelectionStrategy();
		request.genomeSelectionPsmCoverage = params.getGenomeSelectionPsmCoverage();
		request.genomeSelectionStepGuard = params.getGenomeSelectionStepGuard();
		request.genomeSelectionMaxGenomes = params.getGenomeSelectionMaxGenomes();
	}

	static String genomeSelectionSummary(DDAMagOpenSearchParameters params) {
		if (GenomeSelectionCut.isPsmCoverage(params.getGenomeSelectionStrategy())) {
			return "PSM coverage, target " + params.getGenomeSelectionPsmCoverage() + ", step guard "
					+ params.getGenomeSelectionStepGuard() + ", max genomes " + params.getGenomeSelectionMaxGenomes();
		}
		return params.getGenomeSelectionStrategy();
	}

	/**
	 * Round-1 (marker) search settings. Closed by default: a ppm precursor window
	 * with isotope errors 0/1/2 and variable Met oxidation and protein N-terminal
	 * acetylation, the modifications FragPipe's closed marker step uses. Open
	 * mode keeps the shared Da window of {@link #baseSettings}.
	 */
	static MSFraggerHeadlessRunner.OpenSearchSettings round1Settings(DDAMagOpenSearchParameters params,
			File markerTd) {
		MSFraggerHeadlessRunner.OpenSearchSettings settings = baseSettings(params, markerTd);
		// Round 1 only. Overridden here rather than inside baseSettings() so round
		// two, which shares that method, keeps its own outputReportTopN: there the
		// extra candidates feed peak-discovery votes and classifier negatives.
		// Round one has no classifier, so nothing can promote a lower rank, and
		// the reader already discards ranks 2-5 - this stops MSFragger writing
		// them at all.
		settings.outputReportTopN = Math.max(1, params.getRound1OutputReportTopN());
		// Round 1 also writes the Percolator feature table MSBooster consumes.
		// pepXML_pin, not pin: the pepXML this pipeline reads stays exactly as it
		// was, and the pin is additional rather than instead.
		if (params.isStage1RescoringEnabled() && params.isStage1MsboosterEnabled()) {
			settings.outputFormat = "pepXML_pin";
		}
		if (params.isRound1ClosedSearch()) {
			// The selector folds +1/+2 isotope offsets onto the zero family, so
			// isotope-error PSMs are kept rather than lost to the narrow window.
			applyClosedSearch(settings, params.getRound1ClosedPrecursorPpm());
		}
		return settings;
	}

	/**
	 * Closed-stage search of a sample's refined database: the round-1 closed settings
	 * against {@code refined}, {@code closed_stage_topn} candidates per spectrum in every mode.
	 */
	static MSFraggerHeadlessRunner.OpenSearchSettings closedStageSettings(DDAMagOpenSearchParameters params,
			File refined) {
		MSFraggerHeadlessRunner.OpenSearchSettings settings = baseSettings(params, refined);
		settings.outputReportTopN = params.getClosedStageTopN();
		applyClosedSearch(settings, params.getClosedStagePrecursorPpm());
		return settings;
	}

	/** +-ppm precursor window, isotope errors 0/1/2, variable M oxidation and protein N-terminal acetyl. */
	private static void applyClosedSearch(MSFraggerHeadlessRunner.OpenSearchSettings settings, double ppm) {
		settings.precursorMassLower = -ppm;
		settings.precursorMassUpper = ppm;
		settings.precursorMassUnits = 1;
		settings.isotopeError = "0/1/2";
		settings.variableMods = new ArrayList<>(java.util.Arrays.asList("15.9949 M 3", "42.0106 [^ 1"));
	}

	/** Crystal-C refines open-search precursors only; a closed marker search skips it. */
	static boolean round1UsesCrystalC(DDAMagOpenSearchParameters params) {
		return params.isCrystalcEnabled() && !params.isRound1ClosedSearch();
	}

	// --- output writers ----------------------------------------------------

	private static void writeUnmodifiedPsms(UnmodifiedPsmSelector.Result selection, File out) throws IOException {
		writeUnmodifiedPsms(selection, null, null, out);
	}

	/**
	 * Adds one trailing {@code mass_shift} column, empty for zero-family rows and
	 * naming the admitted shift otherwise. Appending rather than inserting keeps
	 * every existing column at its existing index for anything already parsing
	 * this file.
	 */
	private static void writeUnmodifiedPsms(UnmodifiedPsmSelector.Result selection,
			ExpectedPeakSelector.Result expected, Stage1Rescorer.Result rescored, File out) throws IOException {
		out.getParentFile().mkdirs();
		try (PrintWriter writer = new PrintWriter(out, "UTF-8")) {
			writer.println("# calibration_offset_da=" + selection.calibrationOffset + " zero_peak_discovery="
					+ selection.discoveryMethod + " zero_peak_lower_da=" + selection.zeroPeakLowerBoundary
					+ " zero_peak_upper_da=" + selection.zeroPeakUpperBoundary + " zero_peak_local_width_da="
					+ selection.zeroPeakSigma + " excluded_outside_zero_peak=" + selection.excluded
					+ " near_isotope_offsets=" + selection.nearIsotopeOffsets + " fdr_threshold="
					+ selection.fdrThreshold + " fdr_estimator=decoys/targets targets_before_fdr="
					+ selection.targetsBeforeFdr + " decoys_before_fdr=" + selection.decoysBeforeFdr
					+ " targets_after_fdr=" + selection.targetsAfterFdr + " decoys_after_fdr="
					+ selection.decoysAfterFdr + " estimated_fdr=" + selection.estimatedFdr
					+ " expected_shift_targets=" + (expected == null ? 0 : expected.admittedTargets)
					+ " expected_shift_decoys=" + (expected == null ? 0 : expected.admittedDecoys)
					+ " expected_shifts_admitted="
					+ (expected == null ? "" : String.join("|", namesOf(expected)))
					// Which score decided the cut, and on what. A run that fell back
					// used to be indistinguishable from one that rescored.
					+ " stage1_score=" + (rescored == null ? "engine_expect" : "classifier")
					+ " rescoring_features=" + (rescored == null ? 0 : rescored.featureNames.size())
					+ " rescoring_predicted_features=" + (rescored == null ? 0 : rescored.withPredictedFeatures)
					+ " rescoring_targets_after_fdr=" + (rescored == null ? 0 : rescored.targetsAfterFdr)
					+ " rescoring_decoys_after_fdr=" + (rescored == null ? 0 : rescored.decoysAfterFdr));
			writer.println("file\tscan\thit_rank\tpeptide\tprotein\tis_target\tcharge\tmass_diff\tnormalized_mass_diff\t"
					+ "observed_peak_id\tpeak_id\tpeak_isotope_offset\thyperscore\texpect\tq_value\tmass_shift");
			List<OpenSearchPsm> rows;
			if (rescored != null) {
				rows = new ArrayList<>(rescored.accepted);
			} else {
				rows = new ArrayList<>(selection.unmodified);
				if (expected != null) {
					rows.addAll(expected.admitted);
				}
			}
			for (OpenSearchPsm psm : rows) {
				writer.println(psm.getFileName() + "\t" + psm.getScan() + "\t" + psm.getHitRank() + "\t"
						+ psm.getPeptide() + "\t" + psm.getProtein() + "\t" + psm.isTarget() + "\t" + psm.getCharge()
						+ "\t" + psm.getMassDiff() + "\t" + psm.getNormalizedMassDiff() + "\t"
						+ psm.getObservedPeakId() + "\t" + psm.getPeakId() + "\t" + psm.getPeakIsotopeOffset()
						+ "\t" + psm.getHyperscore() + "\t" + psm.getExpect() + "\t" + psm.getQValue()
						+ "\t" + (expected == null ? "" : expected.shiftLabels.getOrDefault(psm, "")));
			}
		}
	}

	private static List<String> namesOf(ExpectedPeakSelector.Result expected) {
		List<String> names = new ArrayList<>();
		for (ExpectedPeakSelector.ExpectedShift shift : expected.admittedShifts()) {
			names.add(shift.name);
		}
		return names;
	}

	/**
	 * The audit trail for stage 1's expected shifts: every UNIMOD shift with PSMs in
	 * its window, admitted or not, with the numbers each criterion was judged on. Written even when
	 * nothing is admitted, because "nothing passed" is the answer that most needs
	 * explaining.
	 */
	private static void writeExpectedPeakVerdicts(ExpectedPeakSelector.Result expected, File out) throws IOException {
		out.getParentFile().mkdirs();
		try (PrintWriter writer = new PrintWriter(out, "UTF-8")) {
			writer.println("# stage-1 expected mass shifts. A shift is admitted only if it is a UNIMOD mass (isotopic labels excluded), has "
					+ "enough target PSMs to test, is target-enriched before any cut, and enough of its accepted "
					+ "peptides were also identified UNMODIFIED in this run. Each admitted shift gets its own "
					+ "q-value cut, computed within that peak and never pooled.");
			writer.println("shift_name\tshift_mass_da\tadmitted\ttargets_in_window\tdecoys_in_window\t"
					+ "window_fdr\taccepted_targets\taccepted_decoys\tpeptide_corroboration\treason");
			for (ExpectedPeakSelector.Verdict v : expected.verdicts) {
				writer.println(v.shift.name + "\t" + v.shift.massDa + "\t" + v.admitted + "\t" + v.targetsInWindow
						+ "\t" + v.decoysInWindow + "\t" + v.windowFdr + "\t" + v.acceptedTargets + "\t"
						+ v.acceptedDecoys + "\t" + v.corroboration + "\t" + v.reason);
			}
		}
	}

	private static void writePeaks(DeltaMassPeakFitter.Result fit, List<OpenSearchPsm> psms, File out)
			throws IOException {
		int[] targetCounts = new int[fit.peaks.size()];
		int[] decoyCounts = new int[fit.peaks.size()];
		Map<Integer, Set<String>> targetFiles = new HashMap<>();
		for (int i = 0; i < fit.observedAssignment.length; i++) {
			int peakId = fit.observedAssignment[i];
			if (peakId < 0) {
				continue;
			}
			OpenSearchPsm psm = psms.get(i);
			if (psm.isTarget()) {
				targetCounts[peakId]++;
				targetFiles.computeIfAbsent(peakId, k -> new HashSet<>()).add(psm.getFileName());
			} else {
				decoyCounts[peakId]++;
			}
		}
		out.getParentFile().mkdirs();
		try (PrintWriter writer = new PrintWriter(out, "UTF-8")) {
			writer.println("peak_id\tcentroid_da\tcanonical_centroid_da\tlower_boundary_da\tupper_boundary_da\t"
					+ "sigma_da\tamplitude\tbackground\tsignal_to_background\tpeak_snr\tr2\tcandidate_count\t"
					+ "discovery_support\ttarget_candidate_count\tdecoy_candidate_count\tfiles_with_target_psms\t"
					+ "persistent_scales\treproducible_files\tdiscovery_method\tisotope_offset\tmono_peak_id");
			for (DeltaMassPeakFitter.Peak peak : fit.peaks) {
				double canonicalCentroid = peak.monoPeakId >= 0 && peak.monoPeakId < fit.peaks.size()
						? fit.peaks.get(peak.monoPeakId).centroid : peak.centroid;
				writer.println(peak.id + "\t" + peak.centroid + "\t" + canonicalCentroid + "\t"
						+ peak.lowerBoundary + "\t" + peak.upperBoundary + "\t" + peak.sigma + "\t" + peak.amplitude + "\t"
						+ peak.background + "\t" + peak.signalToBackground() + "\t" + peak.signalToNoise() + "\t"
						+ peak.r2 + "\t" + peak.psmCount + "\t" + peak.discoverySupport
						+ "\t" + targetCounts[peak.id] + "\t" + decoyCounts[peak.id] + "\t"
						+ targetFiles.getOrDefault(peak.id, Collections.emptySet()).size() + "\t" + peak.persistentScales
						+ "\t" + peak.reproducibleFiles + "\t" + peak.discoveryMethod + "\t" + peak.isotopeOffset
						+ "\t" + peak.monoPeakId);
			}
		}
	}

	/**
	 * {@code psm.getRetentionTime()} is {@code retention_time_sec} straight from
	 * MSFragger's pepXML - seconds, as its own name and
	 * {@code writeVerdicts}'s {@code retention_time_sec} diagnostic column both
	 * say. FlashLFQ's generic identification format expects the opposite:
	 * "Scan Retention Time" in minutes, the same convention MetaMorpheus's own
	 * PSM output uses and the one {@code FlashLfqSettings.toml}'s
	 * {@code MbrRtWindow = 1.5} (1.5 minutes, FlashLFQ's documented default) is
	 * defined in. Feeding it raw seconds does not fail loudly - FlashLFQ still
	 * runs - it just searches for each peptide's elution apex roughly 60x too
	 * far from where it actually eluted, so match-between-runs and even direct
	 * MS/MS-triggered quantification miss almost everywhere. Measured on the
	 * 12mix combined_msfragger_sage run this was still true for: of 225,351
	 * quantified peptide rows, 213,808 (95%) had zero intensity in every
	 * sample, and column 15 of QuantifiedPeptides.tsv was 72% NotDetected and
	 * 25% MSMSIdentifiedButNotQuantified - both exactly what FlashLFQ reports
	 * when it cannot locate the stated retention time within any plausible
	 * chromatographic window.
	 */
	private static double retentionTimeMinutes(OpenSearchPsm psm) {
		return psm.getRetentionTime() / 60.0;
	}

	/**
	 * Writes PSMs in FlashLFQ's own generic identification format (the first 13
	 * columns, exact titles - see {@code FlashLfqPsmReader} and
	 * {@code OpenPsmTsvWriter.pFindFormat}, which this mirrors so the same reader
	 * and the same downstream {@code DdaMagPostAnalysisService} post-analysis work
	 * for this workflow too), plus this workflow's own diagnostic columns appended
	 * after them. Two decisions worth being explicit about:
	 * <ul>
	 * <li><b>Decoys stay in.</b> {@code is_target=false} rows are kept, same as
	 * before - they are the audit trail a genome-level FDR calculation needs.
	 * FlashLFQ itself does not know or care what a decoy is; it only quantifies
	 * whatever sequence+protein it is told about. Downstream biological
	 * reporting ({@link OpenSearchPostAnalysisTsvExporter}) is what excludes them,
	 * the same way it always has.</li>
	 * <li><b>Protein Accession is the full candidate list</b> ({@code getAllProteins()},
	 * {@code ;}-joined), not the single pepXML delegate - matching
	 * {@code PFindPSM.getProtein()} and the reasoning in
	 * {@link OpenSearchPostAnalysisTsvExporter#aggregatePeptides}: razor/protein
	 * inference needs every protein a peptide could belong to, not one arbitrary
	 * pick. The single delegate is kept too, as {@code delegate_protein}, purely
	 * for provenance.</li>
	 * </ul>
	 * The substitution-candidate columns this file used to carry inline
	 * ({@code delta_mass_reassignment}, {@code substitution_verdict},
	 * {@code substitution_candidate}) have moved to {@link #writeSubstitutionForms};
	 * this file now holds exactly one row per PSM, which is what FlashLFQ expects.
	 */
	private static void writeValidatedPsms(List<OpenSearchPsm> psms, Set<OpenSearchPsm> identificationOnly, File out)
			throws IOException {
		out.getParentFile().mkdirs();
		try (PrintWriter writer = new PrintWriter(out, "UTF-8")) {
			writer.println("File Name\tScan number\tBase Sequence\tFull Sequence\tPeptide Monoisotopic Mass\t"
					+ "Scan Retention Time\tPrecursor Charge\tProtein Accession\tScore\tPEP\tPeptide Mass\t"
					+ "Mass difference\tMissed cleavages\t"
					+ "hit_rank\tengine\tengine_support\tlocalized_sequence\tdelegate_protein\tis_target\t"
					+ "normalized_mass_diff\tobserved_peak_id\tpeak_id\tpeak_isotope_offset\thyperscore\texpect\t"
					+ "delta_mass_score_gain\tlocalization_delta_score\tclass_score\tclass_logit\tq_value\t"
					+ "q_value_scope\tvalidation_tier\tstage\tid_group\tshift_state\tshift_site\t"
					+ "shift_site_probability\tshift_explanation\tshift_class");
			for (OpenSearchPsm psm : psms) {
				// Membership of the identification-only result, NOT peakId < 0: a PSM
				// demoted from a rejected peak has a peak id and is still an
				// identification, and labelling it PEAK_VALIDATED would assert a
				// modification the filter explicitly refused.
				boolean unexplained = identificationOnly.contains(psm);
				String fullSequence = OpenSearchPostAnalysisTsvExporter.canonicalModifiedSequenceOf(psm);
				String proteinAccession = joinProteins(psm);
				// "Peptide Monoisotopic Mass" is psm.getFlashLfqMonoisotopicMass(),
				// NOT psm.getCalcNeutralPepMass() - see that method's javadoc. pFind's
				// OpenPsmTsvWriter.addPFindPSM convention (calculated mass here,
				// observed mass in "Peptide Mass") is correct for a closed search,
				// where the calculated mass already includes every applied
				// modification; it silently omits the open delta mass here, sending
				// FlashLFQ to search the wrong region of the MS1 spectrum for every
				// modified peptide (measured: 33% of validated PSMs on the 12mix run).
				// "Peptide Mass" keeps the observed/experimental mass for diagnostics,
				// and the "PEP" slot (this pipeline has no literal posterior error
				// probability) carries the q-value.
				writer.println(psm.getFileName() + "\t" + psm.getScan() + "\t" + psm.getPeptide() + "\t"
						+ fullSequence + "\t" + psm.getFlashLfqMonoisotopicMass() + "\t" + retentionTimeMinutes(psm)
						+ "\t" + psm.getCharge() + "\t" + proteinAccession + "\t"
						+ OpenSearchPostAnalysisTsvExporter.score(psm) + "\t" + psm.getQValue() + "\t"
						+ psm.getPrecursorNeutralMass() + "\t" + psm.getMassDiff() + "\t"
						+ psm.getNumMissedCleavages() + "\t"
						+ psm.getHitRank() + "\t" + psm.getEngine() + "\t" + psm.getEngineSupport() + "\t"
						+ nullSafe(psm.getLocalizedSequence()) + "\t" + nullSafe(psm.getProtein()) + "\t"
						+ psm.isTarget() + "\t" + psm.getNormalizedMassDiff() + "\t" + psm.getObservedPeakId()
						+ "\t" + psm.getPeakId() + "\t" + psm.getPeakIsotopeOffset() + "\t" + psm.getHyperscore()
						+ "\t" + psm.getExpect() + "\t" + psm.getDeltaMassScoreGain() + "\t"
						+ psm.getLocalizationDeltaScore() + "\t" + psm.getClassScore() + "\t" + psm.getClassLogit()
						+ "\t" + psm.getQValue() + "\t" + qValueScope(psm, unexplained)
						+ "\t" + (psm.isClosedStage() ? "CLOSED_VALIDATED"
								: unexplained ? "IDENTIFICATION_ONLY" : "PEAK_VALIDATED")
						+ "\t" + (psm.isClosedStage() ? "closed" : "open") + "\t"
						+ nullSafe(psm.getIdentificationGroup()) + "\t" + shiftColumns(psm));
			}
		}
	}

	/**
	 * {@code shift_state}, {@code shift_site} (1-based), {@code shift_site_probability} and
	 * {@code shift_explanation} (identification v2): this PSM's own fragment-supported localisation,
	 * for localised and accepted-shift PSMs, and for exact-mass PSMs whose isotope offset was an
	 * established shift. "Full Sequence" shows the site of its peptide's representative PSM.
	 * {@code assigned}: the fragments do not place an accepted peak's shift, and the site and explanation
	 * are the UNIMOD modification of the peak's mass on this peptide (no site probability). Empty for
	 * every other PSM. Then {@code shift_class}: {@code modification}, {@code substitution} or
	 * {@code modification_or_substitution} for an explained shift, by what fits its residue
	 * ({@link ShiftLocaliser.Result#explanationClass()}); empty otherwise.
	 */
	static String shiftColumns(OpenSearchPsm psm) {
		ShiftLocaliser.Result loc = psm.getShiftLocalisation();
		if (loc == null || loc.site < 0) {
			return "\t\t\t\t";
		}
		String state = loc.state.name().toLowerCase(java.util.Locale.ROOT);
		int site = psm.getLocalisedShiftSite();
		if (site < 0) {
			return state + "\t\t\t\t";
		}
		return state + "\t" + (site + 1) + "\t"
				+ (Double.isNaN(loc.siteProbability) ? "" : String.format(java.util.Locale.ROOT, "%.4f", loc.siteProbability))
				+ "\t" + loc.explanationText() + "\t" + loc.explanationClassText();
	}

	/** Which walk produced a PSM's q-value. */
	private static String qValueScope(OpenSearchPsm psm, boolean unexplained) {
		if (psm.getIdentificationGroup() != null) {
			// The joint identification walks each evidence group on its own.
			return "JOINT_GROUP_" + psm.getIdentificationGroup();
		}
		if (psm.isClosedStage()) {
			return "CLOSED_SPECTRUM";
		}
		if (unexplained) {
			return psm.getPeakId() < 0 ? "NO_PEAK_GROUP" : "REJECTED_PEAK_GROUP";
		}
		return "LOCAL_PEAK";
	}

	/**
	 * {@code ;}-joined {@code getAllProteins()}, falling back to the single
	 * delegate when a PSM has no alternative-protein list at all - same fallback
	 * {@link OpenSearchPostAnalysisTsvExporter#aggregatePeptides} uses. Per-row,
	 * per-sample, deliberately not reconciled across samples - see the comment at
	 * this method's call site for why that is correct, not an oversight.
	 */
	private static String joinProteins(OpenSearchPsm psm) {
		String[] all = psm.getAllProteins();
		if (all != null && all.length > 0) {
			StringBuilder sb = new StringBuilder();
			for (String protein : all) {
				if (protein == null || protein.isBlank()) {
					continue;
				}
				if (sb.length() > 0) {
					sb.append(";");
				}
				sb.append(protein.trim());
			}
			if (sb.length() > 0) {
				return sb.toString();
			}
		}
		return nullSafe(psm.getProtein());
	}

	/**
	 * Quantifies {@link OpenSearchLayout#substitutionQuantPsms} with FlashLFQ (see
	 * {@link SubstitutionQuantification}): match-between-runs and normalisation off, the peptide and
	 * peak tables kept in {@code sequence_variants/} under their own names, the peptide table annotated
	 * with roles, variant ids, variant/reference ratios and shared main-run features. A missing
	 * executable, an empty PSM file, or a FlashLFQ failure is reported and does not fail the workflow.
	 * When validation did not run, or an older run wrote the PSM file without roles, it is rebuilt from
	 * {@code substitution_forms.tsv}, {@code unexplained_psms.tsv} and {@code validated_psms.tsv}.
	 */
	static File quantifySubstitutionPsms(DDAMagOpenSearchParameters params, File flashLfqExe,
			int flashLfqThreads, File resultRoot, String engineMode, File spectraDir, ProgressCallback cb) {
		try {
			return quantifySubstitutionPsmsOnce(params, flashLfqExe, flashLfqThreads, resultRoot, engineMode,
					spectraDir, cb);
		} catch (Exception e) {
			LOGGER.warn("Substitution quantification failed; the main result is unchanged", e);
			cb.onMessage("WARNING: substitution quantification failed (" + e.getMessage()
					+ "); the main result is unchanged");
			return null;
		}
	}

	private static File quantifySubstitutionPsmsOnce(DDAMagOpenSearchParameters params, File flashLfqExe,
			int flashLfqThreads, File resultRoot, String engineMode, File spectraDir, ProgressCallback cb)
			throws IOException {
		File psmFile = OpenSearchLayout.substitutionQuantPsms(resultRoot, engineMode);
		File quantified = OpenSearchLayout.substitutionQuantifiedPeptides(resultRoot, engineMode);
		File sidePeaks = OpenSearchLayout.substitutionQuantifiedPeaks(resultRoot, engineMode);
		File mainPeaks = OpenSearchLayout.quantifiedPeaks(resultRoot, engineMode);
		if (!SubstitutionQuantification.hasRoles(psmFile)) {
			SubstitutionQuantification.writePsmsFromReports(OpenSearchLayout.substitutionForms(resultRoot, engineMode),
					OpenSearchLayout.unexplainedPsms(resultRoot, engineMode),
					OpenSearchLayout.validatedPsms(resultRoot, engineMode), psmFile);
		}
		int rows = SubstitutionQuantification.variantRows(psmFile);
		if (rows == 0) {
			Files.deleteIfExists(quantified.toPath());
			Files.deleteIfExists(sidePeaks.toPath());
			cb.onMessage("Substitution quantification: no target substitution PSMs outside validated_psms.tsv; "
					+ "FlashLFQ not run");
			return null;
		}
		if (flashLfqExe == null || !flashLfqExe.isFile()) {
			cb.onMessage("Substitution quantification: " + rows + " PSMs in " + psmFile.getName()
					+ "; no FlashLFQ executable supplied, so they were not quantified");
			return null;
		}
		boolean current = mainPeaks.isFile() ? OpenSearchLayout.isSatisfied(quantified, psmFile, mainPeaks)
				: OpenSearchLayout.isSatisfied(quantified, psmFile);
		if (current) {
			cb.onMessage("Substitution quantification: " + quantified.getName() + " is current, skipping FlashLFQ");
			return quantified;
		}
		if (spectraDir == null || !spectraDir.isDirectory()) {
			cb.onMessage("WARNING: substitution quantification skipped; the spectrum folder is not available ("
					+ spectraDir + "). The main result is unchanged");
			return null;
		}
		MetaData metadata = params.getMetadata();
		if (metadata == null || metadata.getRawFileCount() == 0) {
			cb.onMessage("WARNING: substitution quantification skipped; the parameters carry no raw-file metadata. "
					+ "The main result is unchanged");
			return null;
		}
		cb.onMessage("Substitution quantification: FlashLFQ over " + rows + " variant PSMs and their reference PSMs in "
				+ psmFile.getName() + ", match-between-runs and normalisation off");
		SubstitutionQuantification.SpectraFolderGuard guard = SubstitutionQuantification.SpectraFolderGuard
				.take(spectraDir);
		boolean ok;
		try {
			FlashLfqTask task = new FlashLfqTask(flashLfqExe.getAbsolutePath(), flashLfqThreads);
			ok = task.run(metadata, psmFile.getAbsolutePath(), spectraDir.getAbsolutePath(), false, false);
		} finally {
			guard.restore();
		}
		File dir = psmFile.getParentFile();
		File nativePeptides = new File(dir, "QuantifiedPeptides.tsv");
		File nativePeaks = new File(dir, "QuantifiedPeaks.tsv");
		// One UNASSIGNED protein: FlashLFQ's protein table says nothing here.
		Files.deleteIfExists(new File(dir, "QuantifiedProteins.tsv").toPath());
		if (!ok || !nativePeptides.isFile() || nativePeptides.length() == 0) {
			Files.deleteIfExists(nativePeptides.toPath());
			Files.deleteIfExists(nativePeaks.toPath());
			cb.onMessage("WARNING: substitution quantification failed for " + psmFile.getName()
					+ "; the main result is unchanged");
			return null;
		}
		if (nativePeaks.isFile()) {
			Files.move(nativePeaks.toPath(), sidePeaks.toPath(), StandardCopyOption.REPLACE_EXISTING);
		} else {
			Files.deleteIfExists(sidePeaks.toPath());
		}
		int sharing = SubstitutionQuantification.annotate(psmFile, nativePeptides, sidePeaks, mainPeaks);
		Files.move(nativePeptides.toPath(), quantified.toPath(), StandardCopyOption.REPLACE_EXISTING);
		cb.onMessage("Substitution quantification: " + quantified.getName() + " written; " + sharing
				+ " variant peptide(s) share an MS1 feature with a main-run peptide (shares_main_feature)"
				+ (mainPeaks.isFile() ? "" : ", not checked: no main-run QuantifiedPeaks.tsv")
				+ "; final tables are unchanged");
		return quantified;
	}

	/**
	 * One row per PSM that {@link SubstitutionReassigner} flagged (non-empty
	 * {@code substitutionVerdict}), carrying the original call and the candidate
	 * form it was checked against. Joins back to {@code validated_psms.tsv} on
	 * (file, scan, hit_rank) - the same triple that identifies a PSM row there.
	 * Kept separate so a PSM's substitution history never has to be reconciled
	 * with "exactly one row per PSM", which is what FlashLFQ requires of its
	 * identification input.
	 * <p>
	 * {@code variant_id}/{@code variant_class} link a PSM to its row in
	 * {@code sequence_variants.tsv}; empty for decoys, for candidates already in
	 * the searched database, and when variant calling is off.
	 */
	private static void writeSubstitutionForms(List<OpenSearchPsm> psms, SequenceVariantCaller.Result variants,
			File out) throws IOException {
		out.getParentFile().mkdirs();
		try (PrintWriter writer = new PrintWriter(out, "UTF-8")) {
			writer.println("file\tscan\thit_rank\tpeak_id\toriginal_peptide\toriginal_protein\t"
					+ "delta_mass_reassignment\tsubstitution_verdict\tsubstitution_candidate\tvariant_id\t"
					+ "variant_class");
			for (OpenSearchPsm psm : psms) {
				String verdict = psm.getSubstitutionVerdict();
				if (verdict == null || verdict.isBlank()) {
					continue;
				}
				SequenceVariantCaller.Variant variant = variants == null ? null : variants.byPsm.get(psm);
				writer.println(psm.getFileName() + "\t" + psm.getScan() + "\t" + psm.getHitRank() + "\t"
						+ psm.getPeakId() + "\t" + psm.getPeptide() + "\t" + nullSafe(psm.getProtein()) + "\t"
						+ nullSafe(psm.getDeltaMassReassignment()) + "\t" + verdict + "\t"
						+ nullSafe(psm.getSubstitutionCandidate()) + "\t"
						+ (variant == null ? "" : variant.id) + "\t"
						+ (variant == null ? "" : variant.variantClass.name()));
			}
		}
	}

	/**
	 * Complete classifier audit table. Unlike {@code validated_psms.tsv}, this
	 * contains every engine candidate, including non-winning and rejected rows,
	 * so rank pruning, classifier errors and gate decisions can be measured from
	 * one durable artifact rather than reconstructed from pepXML.
	 */
	private static void writeScoredPsms(List<OpenSearchPsm> allCandidates, List<OpenSearchPsm> winners,
			Map<OpenSearchPsm, Double> spectrumQValues, ModificationValidator.Result validation,
			ClosedOpenCascade.ClosedStage closed, File out) throws IOException {
		Set<OpenSearchPsm> winnerSet = Collections.newSetFromMap(new IdentityHashMap<>());
		winnerSet.addAll(winners);
		Set<OpenSearchPsm> validatedSet = Collections.newSetFromMap(new IdentityHashMap<>());
		validatedSet.addAll(validation.validated);
		Set<OpenSearchPsm> peakValidatedSet = Collections.newSetFromMap(new IdentityHashMap<>());
		peakValidatedSet.addAll(validation.peakValidated);
		Set<OpenSearchPsm> unexplainedSet = Collections.newSetFromMap(new IdentityHashMap<>());
		unexplainedSet.addAll(validation.unexplainedIdentifications);
		Map<Integer, Boolean> acceptedPeaks = new HashMap<>();
		for (ModificationValidator.PeakVerdict verdict : validation.verdicts) {
			acceptedPeaks.put(verdict.peakId, verdict.accepted);
		}

		out.getParentFile().mkdirs();
		try (PrintWriter writer = new PrintWriter(out, "UTF-8")) {
			// engine first after the identifiers: without it a mixed-engine run
			// cannot be audited, and "did Sage add anything" is unanswerable.
			writer.println("file\tscan\thit_rank\tengine\tengine_support\tpeptide\tlocalized_sequence\tprotein\tis_target\tcharge\t"
					+ "precursor_neutral_mass\tretention_time_sec\tmass_diff\tnormalized_mass_diff\tobserved_peak_id\tpeak_id\tpeak_centroid\t"
					+ "peak_isotope_offset\tdistance_to_peak_centroid\tdistance_in_sigmas\tlocal_density\t"
					+ "hyperscore\tnextscore\texpect\tscore_with_delta_mass\tscore_without_delta_mass\t"
					+ "delta_mass_score_gain\tsecond_best_score_with_delta_mass\tlocalization_delta_score\t"
					+ "ions_with_delta_mass\tions_without_delta_mass\tlocalization_site_count\tptm_mass\t"
					+ "best_positions\tnumber_best_positions\tposition_scores\t"
					+ "class_score\tclass_logit\tis_spectrum_winner\tspectrum_q_value\t"
					+ "peak_hypothesis_accepted\tlocal_q_value\tpeak_validated\t"
					+ "identification_only\tvalidated\tstage\tclosed_accepted");
			List<OpenSearchPsm> rows = new ArrayList<>(allCandidates);
			if (closed != null) {
				rows.addAll(closed.winners);
				for (OpenSearchPsm psm : closed.winners) {
					winnerSet.add(psm);
					spectrumQValues.put(psm, psm.getQValue());
				}
			}
			for (OpenSearchPsm psm : rows) {
				boolean winner = winnerSet.contains(psm);
				boolean peakHypothesisAccepted = !psm.isClosedStage()
						&& acceptedPeaks.getOrDefault(psm.getPeakId(), false);
				boolean closedAccepted = closed != null
						&& closed.acceptedBySpectrum.containsKey(ClosedOpenCascade.spectrumKey(psm));
				// Every spectrum winner goes through some FDR walk - the per-peak one if
				// its peak was accepted, otherwise the identification walk for rejected-
				// peak PSMs or the separate one for PSMs outside every peak - so a winner
				// always carries a q-value. With the closed stage, an open winner of a
				// spectrum the closed stage accepted goes through none.
				boolean hasLocalQ = winner && (psm.isClosedStage() || !closedAccepted);
				writer.println(psm.getFileName() + "\t" + psm.getScan() + "\t" + psm.getHitRank() + "\t"
						+ psm.getEngine() + "\t" + psm.getEngineSupport() + "\t"
						+ psm.getPeptide() + "\t" + nullSafe(psm.getLocalizedSequence()) + "\t" + psm.getProtein()
						+ "\t" + psm.isTarget() + "\t" + psm.getCharge() + "\t" + psm.getPrecursorNeutralMass()
						+ "\t" + psm.getRetentionTime() + "\t" + psm.getMassDiff() + "\t" + psm.getNormalizedMassDiff()
						+ "\t" + psm.getObservedPeakId() + "\t" + psm.getPeakId()
						+ "\t" + psm.getPeakCentroid() + "\t" + psm.getPeakIsotopeOffset() + "\t"
						+ psm.getDistanceToPeakCentroid() + "\t" + psm.getDistanceInSigmas() + "\t"
						+ psm.getLocalDensity() + "\t" + psm.getHyperscore() + "\t" + psm.getNextscore() + "\t"
						+ psm.getExpect() + "\t" + psm.getScoreWithDeltaMass() + "\t"
						+ psm.getScoreWithoutDeltaMass() + "\t" + psm.getDeltaMassScoreGain() + "\t"
						+ psm.getSecondBestScoreWithDeltaMass() + "\t" + psm.getLocalizationDeltaScore() + "\t"
						+ psm.getIonsWithDeltaMass() + "\t" + psm.getIonsWithoutDeltaMass() + "\t"
						+ psm.getLocalizationSiteCount() + "\t" + ptmAttribute(psm, "ptm_mass") + "\t"
						+ ptmAttribute(psm, "best_positions") + "\t" + ptmAttribute(psm, "number_best_positions")
						+ "\t" + ptmAttribute(psm, "position_scores") + "\t" + psm.getClassScore() + "\t"
						+ psm.getClassLogit()
						+ "\t" + winner + "\t" + valueOrBlank(spectrumQValues.get(psm)) + "\t"
						+ peakHypothesisAccepted + "\t" + (hasLocalQ ? valueOrBlank(psm.getQValue()) : "") + "\t"
						+ peakValidatedSet.contains(psm) + "\t" + unexplainedSet.contains(psm) + "\t"
						+ validatedSet.contains(psm) + "\t" + (psm.isClosedStage() ? "closed" : "open") + "\t"
						+ closedAccepted);
			}
		}
	}

	private static void writeVerdicts(ModificationValidator.Result validation, File out) throws IOException {
		out.getParentFile().mkdirs();
		try (PrintWriter writer = new PrintWriter(out, "UTF-8")) {
			writer.println("# mass_class: ZERO is the peak nearest 0 Da and is always accepted; ESTABLISHED passed the");
			writer.println("# evidence rule on its rank-1 PSMs of 13+ residues: 50 or more of them, at most 5% decoys,");
			writer.println("# and a long-peptide share of at least 0.75 of that of all other non-zero peaks;");
			writer.println("# UNEXPLAINED failed it - its PSMs are identifications with no modification claim.");
			writer.println("# A peak that passes but whose localised PSMs mostly sit where a substitution fits is");
			writer.println("# withdrawn (accepted=false, see reason and peak_site_shares.tsv).");
			writer.println("# fragment_supported is recorded from the delta-mass evidence gate, not used to accept.");
			writer.println("# Decoy fractions are screening statistics, not calibrated: peaks are defined by fitting");
			writer.println("# a target-dominated histogram. FDR control comes from the walks using 1:1 decoys.");
			writer.println("# peak_id=-1 is deliberately absent: it is a mixture, not a modification hypothesis.");
			writer.println("# accepted=false does NOT discard that peak's PSMs. They are cut by an identification");
			writer.println("# FDR walk over rejected-peak PSMs, separate from the walk over peak_id=-1;");
			writer.println("# psms_retained counts peak-validated PSMs only, so it is 0 for a rejected peak.");
			writer.println("# identification_only_retained=" + validation.unexplainedIdentifications.size()
					+ " (PSMs with no accepted mass explanation, retained by identification FDR)");
			writer.println("peak_id\tcentroid_da\tmass_class\ttargets\tdecoys\tdecoy_fraction\t"
					+ "decoy_upper_bound\tbackground_decoy_fraction\tfragment_supported\taccepted\treason\t"
					+ "psms_retained\tlong_peptide_psms\tlong_peptide_decoy_fraction\tlong_peptide_share_ratio");
			for (ModificationValidator.PeakVerdict verdict : validation.verdicts) {
				writer.println(verdict.peakId + "\t" + verdict.centroid + "\t" + verdict.massClass + "\t"
						+ verdict.targets + "\t" + verdict.decoys + "\t"
						+ verdict.decoyFraction + "\t" + verdict.decoyUpperBound + "\t"
						+ verdict.backgroundDecoyFraction + "\t" + verdict.fragmentSupported + "\t"
						+ verdict.accepted + "\t" + verdict.reason + "\t" + verdict.psmsRetained + "\t"
						+ verdict.longPeptidePsms + "\t" + verdict.longPeptideDecoyFraction + "\t"
						+ verdict.longPeptideShareRatio);
			}
		}
	}

	/**
	 * Tolerance for matching a peak centroid to a candidate modification or
	 * substitution mass. Matches {@link SubstitutionReassigner.Config#massToleranceDa}'s
	 * default and its stated rationale: fitted peaks in this workflow have sigma
	 * around 2 mDa, so 5 mDa admits the peak and little else.
	 */
	private static final double POTENTIAL_MODIFICATION_TOLERANCE_DA = 0.005;

	/**
	 * One row per <em>accepted</em> delta-mass peak: what its mass might be, not
	 * just whether it passed FDR.
	 * <p>
	 * Two independent lookups run against the same centroid - {@link
	 * ModificationMassDatabase}, the bundled 586-entry chemical-modification
	 * reference (Kai Cheng's direction, 2026-09-22: "Use the bundled 586-entry set
	 * in current stage"), and {@link SubstitutionMassTable}, the single-residue
	 * substitution table {@link SubstitutionReassigner} already uses. Neither
	 * lookup is gated on peak-level {@code substitutionShaped}: a peak can have
	 * <em>some</em> PSMs whose localization is substitution-consistent without
	 * clearing that per-peak threshold, and {@code candidate_substitutions} reports
	 * whatever share was found, via {@code substitutions.peaks} - which, unlike
	 * {@code substitutions.substitutionPeakIds}, is not restricted to shaped peaks.
	 * <p>
	 * A peak that <em>did</em> clear the substitution-shaped threshold is withdrawn
	 * from {@code gate.acceptedPeakIds} before {@code validation} is computed (see
	 * this method's call site), so it is never {@code accepted} here and never
	 * appears in this file - it is already reported as a substitution, not a
	 * modification, in {@code substitution_peaks.tsv}.
	 */
	private static void writePotentialModifications(DeltaMassPeakFitter.Result fit,
			SubstitutionReassigner.Result substitutions, ModificationValidator.Result validation, File out)
			throws IOException {
		Map<Integer, DeltaMassPeakFitter.Peak> peaksById = new HashMap<>();
		for (DeltaMassPeakFitter.Peak peak : fit.peaks) {
			peaksById.put(peak.id, peak);
		}
		Map<Integer, SubstitutionReassigner.PeakSummary> substitutionsById = new HashMap<>();
		for (SubstitutionReassigner.PeakSummary s : substitutions.peaks) {
			substitutionsById.put(s.peakId, s);
		}

		out.getParentFile().mkdirs();
		try (PrintWriter writer = new PrintWriter(out, "UTF-8")) {
			writer.println("# One row per accepted delta-mass peak (see modification_verdicts.tsv for every peak,");
			writer.println("# accepted or not). candidate_modifications lists every entry in the bundled 586-entry");
			writer.println("# modifications.xml reference within " + POTENTIAL_MODIFICATION_TOLERANCE_DA
					+ " Da of the peak's fitted centroid, nearest first, as name[mass_error_da].");
			writer.println("# candidate_substitutions lists every single amino-acid substitution within the same");
			writer.println("# tolerance, plus the share of this peak's PSMs whose localization is consistent with a");
			writer.println("# substitution and the share of those already in the searched database - see");
			writer.println("# SubstitutionReassigner. A peak judged substitution-shaped there is withdrawn before");
			writer.println("# acceptance and never appears here; it is reported instead in substitution_peaks.tsv.");
			writer.println("# best_candidate_* picks the smaller-|error| of the two lookups' nearest match; UNEXPLAINED");
			writer.println("# means neither lookup found anything within tolerance.");
			writer.println("peak_id\tcentroid_da\tsigma_da\ttarget_candidate_count\tdecoy_candidate_count\t"
					+ "decoy_fraction\tdecoy_upper_bound\taccepted\tpeak_snr\tr2\tdiscovery_support\t"
					+ "candidate_modifications\tcandidate_substitutions\tbest_candidate_type\tbest_candidate_name\t"
					+ "mass_error_da");
			for (ModificationValidator.PeakVerdict verdict : validation.verdicts) {
				if (!verdict.accepted) {
					continue;
				}
				DeltaMassPeakFitter.Peak peak = peaksById.get(verdict.peakId);
				if (peak == null) {
					continue;
				}
				List<ModificationMassDatabase.Match> modMatches = ModificationMassDatabase.matching(peak.centroid,
						POTENTIAL_MODIFICATION_TOLERANCE_DA);
				List<SubstitutionMassTable.Substitution> subMatches = SubstitutionMassTable.matching(peak.centroid,
						POTENTIAL_MODIFICATION_TOLERANCE_DA);

				String bestType = "UNEXPLAINED";
				String bestName = "";
				String bestError = "";
				double bestAbsError = Double.POSITIVE_INFINITY;
				if (!modMatches.isEmpty()) {
					bestType = "MODIFICATION";
					bestName = modMatches.get(0).entry.title;
					bestAbsError = Math.abs(modMatches.get(0).massErrorDa);
					bestError = String.format(java.util.Locale.ROOT, "%.5f", modMatches.get(0).massErrorDa);
				}
				if (!subMatches.isEmpty()) {
					double subError = peak.centroid - subMatches.get(0).massDelta;
					if (Math.abs(subError) < bestAbsError) {
						bestType = "SUBSTITUTION";
						bestName = subMatches.get(0).toString();
						bestError = String.format(java.util.Locale.ROOT, "%.5f", subError);
					}
				}

				writer.println(verdict.peakId + "\t" + peak.centroid + "\t" + peak.sigma + "\t" + verdict.targets
						+ "\t" + verdict.decoys + "\t" + verdict.decoyFraction + "\t" + verdict.decoyUpperBound + "\t"
						+ verdict.accepted + "\t" + peak.signalToNoise() + "\t" + peak.r2 + "\t"
						+ peak.discoverySupport + "\t" + formatModifications(modMatches)
						+ "\t" + formatSubstitutions(subMatches, substitutionsById.get(verdict.peakId)) + "\t"
						+ bestType + "\t" + bestName + "\t" + bestError);
			}
		}
	}

	/** Nearest five, {@code title[+/-mass_error_da]}, {@code ;}-joined. */
	private static String formatModifications(List<ModificationMassDatabase.Match> matches) {
		StringBuilder sb = new StringBuilder();
		int shown = 0;
		for (ModificationMassDatabase.Match m : matches) {
			if (shown >= 5) {
				break;
			}
			if (sb.length() > 0) {
				sb.append(";");
			}
			sb.append(m.entry.title).append('[')
					.append(String.format(java.util.Locale.ROOT, "%+.5f", m.massErrorDa)).append(']');
			shown++;
		}
		return sb.toString();
	}

	/** Every substitution within tolerance, {@code ,}-joined, plus the peak's localization/DB support if known. */
	private static String formatSubstitutions(List<SubstitutionMassTable.Substitution> matches,
			SubstitutionReassigner.PeakSummary summary) {
		if (matches.isEmpty()) {
			return "";
		}
		StringBuilder sb = new StringBuilder();
		for (SubstitutionMassTable.Substitution s : matches) {
			if (sb.length() > 0) {
				sb.append(",");
			}
			sb.append(s.toString());
		}
		if (summary != null) {
			sb.append(String.format(java.util.Locale.ROOT, " (%.1f%% localization-consistent, %.1f%% of those in DB)",
					100 * summary.substitutionFraction(), 100 * summary.inDatabaseFraction()));
		}
		return sb.toString();
	}

	/** Writes a combined candidate's score onto one of its engine rows, if that engine proposed it. */
	private static void applyScore(OpenSearchPsm psm, double score, double logit, String engineSupport) {
		if (psm != null) {
			psm.setClassScore(score);
			psm.setClassLogit(logit);
			psm.setEngineSupport(engineSupport);
		}
	}

	private static String nullSafe(String value) {
		return value == null ? "" : value;
	}

	private static String valueOrBlank(Double value) {
		return value == null || Double.isNaN(value) ? "" : value.toString();
	}

	private static String ptmAttribute(OpenSearchPsm psm, String name) {
		return nullSafe(psm.getPtmResultAttributes().get(name));
	}

	private static int countLines(File file) {
		try {
			return (int) Files.lines(file.toPath()).count();
		} catch (IOException e) {
			return 0;
		}
	}

	/**
	 * FragPipe MAG equivalent of {@code cleanupMarkerSearchIntermediates}: after
	 * every search has written its pepXML, drop MSFragger {@code .pepindex} files
	 * under {@code search_results/} and any fallback marker index beside
	 * {@code marker_target_decoy.fasta} in the engine folder. They sit beside the
	 * FASTA MSFragger searched; resume still skips from the pepXML. Paths are
	 * resolved without {@link OpenSearchLayout#searchResultsRoot} so this step
	 * does not migrate folders. The shared MAG-database {@code Marker4_td.fasta}
	 * index is outside this tree and is left in place.
	 */
	static void cleanupDatabaseIndexes(File resultRoot, List<String[]> experiments, ProgressCallback cb) {
		if (resultRoot == null || cb == null) {
			return;
		}
		deletePepindexFiles(resultRoot, cb);
		deletePepindexTree(new File(resultRoot, OpenSearchLayout.SEARCH_RESULTS_DIR), cb);
	}

	private static void deletePepindexTree(File directory, ProgressCallback cb) {
		if (directory == null || !directory.isDirectory()) {
			return;
		}
		deletePepindexFiles(directory, cb);
		File[] children = directory.listFiles(File::isDirectory);
		if (children == null) {
			return;
		}
		for (File child : children) {
			deletePepindexTree(child, cb);
		}
	}

	private static void deletePepindexFiles(File directory, ProgressCallback cb) {
		if (directory == null || !directory.isDirectory()) {
			return;
		}
		File[] files = directory.listFiles(File::isFile);
		if (files == null) {
			return;
		}
		for (File file : files) {
			if (!file.getName().toLowerCase(java.util.Locale.ROOT).endsWith(".pepindex")) {
				continue;
			}
			try {
				if (Files.deleteIfExists(file.toPath())) {
					cb.onMessage("Deleted MSFragger database index: " + file.getAbsolutePath());
				}
			} catch (IOException e) {
				cb.onMessage("WARNING: could not delete MSFragger database index " + file.getAbsolutePath()
						+ ": " + e.getMessage());
			}
		}
	}

	private static String baseName(String path) {
		String name = new File(path).getName();
		int dot = name.lastIndexOf('.');
		return dot > 0 ? name.substring(0, dot) : name;
	}

	/**
	 * The folder the spectra came from, used by the identification-summary
	 * writer to report per-file scan counts. Null when it cannot be determined,
	 * which the writer tolerates.
	 */
	/**
	 * The folder the spectrum-reading tools (MSBooster, Crystal-C, FlashLFQ) are pointed at: the
	 * inputs' folder for mzML inputs, {@code spectra_mzml/} for .raw and .d inputs
	 * ({@link OpenSearchSpectra}).
	 */
	static File readableSpectraDir(List<String[]> experiments, File resultRoot) throws IOException {
		if (experiments == null || experiments.isEmpty()) {
			return null;
		}
		List<File> inputs = new ArrayList<>(experiments.size());
		for (String[] experiment : experiments) {
			inputs.add(new File(experiment[0]));
		}
		return OpenSearchSpectra.spectraDir(inputs, resultRoot);
	}

	private static File spectraDirOf(List<String[]> experiments) {
		if (experiments == null || experiments.isEmpty()) {
			return null;
		}
		File first = new File(experiments.get(0)[0]);
		return first.getParentFile();
	}

	/** Under the result root: {@code ms2Count.txt} and {@code rt.txt} written from the input mzML files. */
	static final String SPECTRA_SUMMARY_DIR = MagResultLayoutUtil.SPECTRA_SUMMARY_DIR_NAME;

	/**
	 * The folder holding {@code ms2Count.txt} (MS2 scans per run) for {@code final_summary.tsv}'s
	 * identification rate. pFind writes it in the result's {@code spectra/} folder. A run started
	 * from mzML converts nothing, so the same two files are written there rather than beside the
	 * inputs. An existing count, including one left in {@code spectra_summary/}, is used as it is.
	 * Falls back to the input folder when the inputs are not all mzML or cannot be read.
	 */
	static File spectraSummaryDir(List<String[]> experiments, File resultRoot, ProgressCallback cb) {
		File inputDir = spectraDirOf(experiments);
		if (MagResultLayoutUtil.hasMs2Count(inputDir)) {
			return inputDir;
		}
		if (experiments == null || experiments.isEmpty() || resultRoot == null) {
			return inputDir;
		}
		String[] files = new String[experiments.size()];
		for (int i = 0; i < files.length; i++) {
			files[i] = OpenSearchSpectra.readable(new File(experiments.get(i)[0]), resultRoot).getPath();
		}
		if (!MzmlSummaryFileWriter.allMzml(files)) {
			cb.onMessage("Spectra summary: inputs are not all mzML; final_summary.tsv uses identified PSMs as the "
					+ "MS/MS count");
			return inputDir;
		}
		File dir = MagResultLayoutUtil.spectraCountDir(resultRoot);
		if (MagResultLayoutUtil.hasMs2Count(dir)) {
			cb.onMessage("Spectra summary: using " + dir.getAbsolutePath());
			return dir;
		}
		try {
			if (MzmlSummaryFileWriter.ensure(dir, files)) {
				cb.onMessage("Spectra summary: ms2Count.txt and rt.txt in " + dir.getAbsolutePath());
				return dir;
			}
			// A partial file would be read as complete next time.
			new File(dir, "ms2Count.txt").delete();
			new File(dir, "rt.txt").delete();
			cb.onMessage("WARNING: Spectra summary: an mzML file is missing or empty; final_summary.tsv uses "
					+ "identified PSMs as the MS/MS count");
		} catch (IOException e) {
			new File(dir, "ms2Count.txt").delete();
			new File(dir, "rt.txt").delete();
			cb.onMessage("WARNING: Spectra summary could not be written (" + e.getMessage()
					+ "); final_summary.tsv uses identified PSMs as the MS/MS count");
		}
		return inputDir;
	}
}
