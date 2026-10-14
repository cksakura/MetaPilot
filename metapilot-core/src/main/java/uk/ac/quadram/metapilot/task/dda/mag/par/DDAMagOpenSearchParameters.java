package uk.ac.quadram.metapilot.task.dda.mag.par;

import uk.ac.quadram.metapilot.core.msfragger.MSFraggerHeadlessRunner;
import java.util.ArrayList;
import java.util.List;

import uk.ac.quadram.metapilot.core.opensearch.DeltaMassPeakFitter;
import uk.ac.quadram.metapilot.core.opensearch.ExpectedPeakSelector;
import uk.ac.quadram.metapilot.core.opensearch.Stage1Rescorer;
import uk.ac.quadram.metapilot.core.opensearch.UnmodifiedPsmSelector;
import uk.ac.quadram.metapilot.core.opensearch.ModificationValidator;
import uk.ac.quadram.metapilot.task.mag.MagDbItem;
import uk.ac.quadram.metapilot.task.par.MetaBaseParameter;
import uk.ac.quadram.metapilot.task.par.MetaConstants;
import uk.ac.quadram.metapilot.task.par.MetaData;

/**
 * Parameters for the standalone DDA MAG open-search workflow.
 * <p>
 * The parameter file is the single source of truth for GUI, CLI and HPC runs
 * alike, matching the convention of the other MAG workflows — there are no
 * CLI-only flags.
 * <p>
 * Note there is no "search mode" setting: <em>both</em> rounds are open
 * searches by construction. Round 1 uses an open search but consumes only its
 * unmodified (zero delta-mass peak) PSMs for genome discovery, because an open
 * search diverts modified spectra to their true peptide with a mass shift
 * rather than force-matching them to a wrong peptide inside a narrow window,
 * leaving a purer set to map genomes from.
 */
public class DDAMagOpenSearchParameters extends MetaBaseParameter {

	/** ExperimentData format: [filePath, experimentName, fraction]. */
	private List<String[]> experimentData = new ArrayList<>();

	private MagDbItem magDbItem;
	private String magDbCatalogueId = "";
	private String magDbIdentifier = "MGYG";
	private String magDbVersion = "";
	private String genomeSelectorMode = MagDbItem.MARKER2;
	private String genomeSelectionStrategy = "Focus";

	/** MSFragger alone. The only mode currently implemented. */
	public static final String ENGINES_MSFRAGGER = "MSFragger only";
	/**
	 * MSFragger + Sage, pooled for union coverage with engine agreement as a
	 * classifier feature. Present so the benefit of a second engine can be
	 * measured against the single-engine baseline rather than assumed — on a
	 * refined, genome-specific database the gain may not justify the cost, and
	 * that is an empirical question. Not yet implemented.
	 */
	public static final String ENGINES_MSFRAGGER_AND_SAGE = "MSFragger + Sage";

	private String searchEngines = ENGINES_MSFRAGGER;

	/** Samples searched concurrently; 0 lets the runner choose from the core count. */
	private int parallelTasks = 0;
	/** -Xmx for the MSFragger JVM. MSFragger is memory hungry on open searches. */
	private String msfraggerMaxHeap = "32g";

	// --- open-search window ------------------------------------------------
	private double precursorMassLower = -150.0;
	private double precursorMassUpper = 500.0;
	/**
	 * PSM candidates MSFragger reports per spectrum. MSFragger's own default is
	 * 1 (top hit only), which would defeat this pipeline's requirement that the
	 * classifier, not the engine, decides between candidates.
	 */
	private int outputReportTopN = 5;
	/**
	 * Candidates reported per spectrum in ROUND 1 only. One, deliberately.
	 * <p>
	 * Round 1 has no classifier - it cannot have one without circularity, since a
	 * well-trained classifier wants the refined database that round 1 exists to
	 * produce - so nothing can ever promote a lower-ranked candidate over the
	 * engine's own ordering, and {@code OpenSearchPepxmlReader.setMaxHitRank(1)}
	 * in the runner already discards ranks 2-5 on read. This setting stops
	 * MSFragger writing them at all.
	 * <p>
	 * Measured on nine 12mix runs: MSFragger's argmax(hyperscore) is rank 1 for
	 * 1,256,620 of 1,256,636 spectra (99.999%); at rank 2 only 8.9% of target hits
	 * map to a true genome against 23.3% at rank 1, falling to 3.1% by rank 5; and
	 * even a perfect oracle rescorer allowed to promote any of ranks 2-5 could add
	 * at most 4.95%. The cost of writing them was 14 GB of round-1 pepXML for nine
	 * runs instead of ~2.6 GB.
	 * <p>
	 * Round 2 keeps {@link #outputReportTopN} unchanged: there the extra
	 * candidates feed spectrum-balanced peak-discovery votes and the
	 * semi-supervised classifier's negatives, which this measurement does not test.
	 */
	// REVIEWED 2026-09-17 (Claude session, source-verified against the javadoc above): rank
	// policy re-confirmed, no open question. See Claude outputs/opensearch_optimization_review.md.
	private int round1OutputReportTopN = 1;

	public static final String ROUND1_SEARCH_CLOSED = "closed";
	public static final String ROUND1_SEARCH_OPEN = "open";
	/**
	 * How round 1 searches the marker database: {@code closed} (default) or {@code open}.
	 * <p>
	 * Closed: +-{@link #round1ClosedPrecursorPpm} ppm, {@code isotope_error = 0/1/2}, variable M
	 * oxidation and protein N-terminal acetyl, as the FragPipe closed workflow; Crystal-C is not run
	 * on it. Open: the round-2 window and {@code isotope_error = 0}, with Crystal-C.
	 * <p>
	 * Round 1 keeps only the zero family plus allow-list shifts, which a closed search delivers. On
	 * 12mix (9 runs) and PXD050781 (21 runs) the per-run genome lists of a closed and an open marker
	 * search with the same selector carry the same spiked genomes (misses 6 and 5 of 108, all but one
	 * the marker-poor <i>E. casseliflavus</i>; PXD050781 target in every run) from the same marker
	 * peptides (6,215 and 6,252 per run), and the closed MSFragger search takes about a fifth of the
	 * time. The closed lists held more non-spiked genomes over all runs (32 against 20 on 12mix), but
	 * they came from FragPipe's processing (Percolator, picked protein filter), not this rescoring.
	 */
	private String round1SearchMode = ROUND1_SEARCH_CLOSED;
	/** Round-1 closed precursor window, +- ppm. */
	private double round1ClosedPrecursorPpm = 20.0;

	/**
	 * Closed search first, open search on the remaining spectra
	 * ({@code docs/developer/closed-then-open-identification-plan.md}). Each sample's refined
	 * database is also searched closed (+-{@link #closedStagePrecursorPpm} ppm, isotope errors
	 * 0/1/2, variable M oxidation and protein N-terminal acetyl); closed PSMs accepted at 1% FDR
	 * are identifications, and the open search identifies only the spectra they leave, as
	 * unmodified or known-shift PSMs. Unexplained shifts are reported but are not
	 * identifications, the fragment-evidence gate no longer folds rejected peaks onto zero, and
	 * allow-list masses are not withdrawn as substitutions. Off restores the open-only workflow.
	 */
	private boolean closedStageEnabled = true;
	/** Closed-stage precursor window, +- ppm. */
	private double closedStagePrecursorPpm = 20.0;

	// --- stage 1: unmodified selection -------------------------------------
	/** Zero-peak membership half-width, in fitted sigmas. */
	private double zeroPeakSigmas = 3.0;
	/** Half-width of the window searched for the zero peak (Da). */
	private double zeroPeakSearchWindow = 0.2;
	/**
	 * Absolute lower bound on peak membership half-width (Da).
	 * <p>
	 * Delta-mass peaks have a sharp core and heavy tails, so a pure 3-sigma rule
	 * cut at ~3 mDa discards real members. The relevant scale is instrument mass
	 * accuracy, not the fitted core width.
	 */
	private double membershipFloorDa = 0.010;

	// --- delta-mass peak fitting -------------------------------------------
	/**
	 * Minimum PSMs in a segment before a Gaussian is fitted to it.
	 * <p>
	 * Defaulted from the binomial gate's own detection floor: with 1:1
	 * target:decoy a peak of n PSMs has best-case p = 0.5^n, so a peak below this
	 * size cannot be reported however clean it is. Fitting smaller segments both
	 * wastes time and inflates the multiple-testing burden on real peaks.
	 */
	private int peakMinArea = 15;
	/** Hard cap on Gaussian fits per run; largest segments kept. 0 = no cap. */
	private int maxPeaksToFit = 500;
	/** Iteration ceiling for one Gaussian fit; non-peaks otherwise run to the default limit. */
	private int maxFitIterations = 200;
	/** Minimum raw bins required to support a four-parameter peak fit. */
	private int peakMinFitBins = 7;
	/** Smallest identifiable fitted sigma, expressed in histogram-bin units. */
	private double peakMinSigmaBins = 0.5;
	/** Minimum amplitude/sqrt(background+1) for a fitted peak. */
	private double peakMinSnr = 5.0;
	/** Absolute minimum separation between distinct fits, in histogram bins. */
	private double peakMinSeparationBins = 2.0;
	/** Width-scaled minimum separation, in combined fitted sigmas. */
	private double peakMinSeparationSigmas = 3.0;

	/** 0 disables; see DeltaMassPeakFitter.Config.kdeMaxSigmaMultiple. */
	private double kdeMaxSigmaMultiple = 0.0;

	/** 0 disables; see DeltaMassPeakFitter.Config.kdeMinProminenceRatio. */
	private double kdeMinProminenceRatio = 0.0;
	/** KDE is the primary discovery model; Gaussian fitting remains a compatibility option. */
	private boolean useKdePeakDiscovery = true;
	private double kdeBandwidthDa = 0.001;
	private String kdeScaleMultipliers = "1/2/4";
	private int kdeMinPersistentScales = 2;
	private double kdeMaxMembershipDa = 0.050;
	private int kdeMinReproducibleFiles = 2;
	private double kdeMinFileSupport = 3.0;
	private int kdeMaxModes = 2000;

	/**
	 * Seed the stage-2 classifier's confident positives only from PSMs inside a
	 * peak the filter admitted. Scoring still covers every candidate.
	 * <p>
	 * <b>Default off, and it should stay off unless you are testing it.</b> When
	 * on, every training positive is by construction inside a validated peak, so
	 * the model never sees a confident PSM that has no accepted mass explanation
	 * and has no basis on which to rank one. On the 2026-09-13 MSFragger-only run
	 * that produced a degenerate classifier:
	 * <pre>
	 *   group           winners   decoys          class_logit      within-group AUC
	 *   accepted_peak   850,268    4.23%    1.110 .. 12.001            0.9733
	 *   rejected_peak    30,126   37.09%    0.318 ..  4.961            0.5035
	 *   no_peak         373,267   35.94%    0.318 ..  1.110            0.4986
	 * </pre>
	 * The score ranges do not overlap: everything in a validated peak scores above
	 * 1.110 and everything else below, with no discrimination whatever outside -
	 * two AUCs indistinguishable from a coin flip. The identification-only walk
	 * then retained 11 PSMs out of 403,393.
	 * <p>
	 * The mechanism is label leakage: {@code in_validated_peak} selected the
	 * training positives <em>and</em> was an input feature, so the model separated
	 * its training set by reading that one column. It was invisible while the
	 * binomial gate accepted nearly every peak, because the feature then carried
	 * almost no label information. A filter that genuinely rejects peaks makes the
	 * feature align with the label, and the model collapses onto it.
	 * {@link uk.ac.quadram.metapilot.core.opensearch.OpenSearchFeatures#withoutLabelDefiningFeatures}
	 * now removes the feature whenever this is on, but the deeper bias - positives
	 * drawn only from inside peaks - remains, which is why the default is off.
	 */
	private boolean restrictPositivesToValidatedPeaks = false;
	/**
	 * Score one instance per (spectrum, peptide) carrying both engines' evidence,
	 * instead of one instance per engine row calibrated onto a common scale.
	 * <p>
	 * Only meaningful with two engines; a MSFragger-only run takes the same path
	 * either way. Default on: measured on the 2026-09-10 combined run the
	 * calibration it replaces is not achieving what it was designed to achieve -
	 * at equal calibrated score Sage rows are about 1.8x more likely to be decoys
	 * (9.47% vs 5.17% in the 1-3 band), because Sage's classifier sees 25 features
	 * against MSFragger's 33 and its engine-local q-value is correspondingly
	 * optimistic. Kept as a switch so the two can be compared on one dataset.
	 */
	private boolean combinedPsmInstances = true;

	/**
	 * Round-two candidates per spectrum. Five.
	 * <p>
	 * This was 1, on an earlier 12mix measurement recorded as "at every matched
	 * entrapment budget rank-1 alone recovered more true PSMs (-1.0% to -6.8% for
	 * top-5)". Re-measured on 2026-09-12 against the 12mix genome ground truth,
	 * that does not reproduce. At a matched 0.49% decoy budget:
	 *
	 * <pre>
	 *   rank-1 only   686,076 targets   80,958 peptides   85.5% purity
	 *   ranks 1-5     791,078 targets   85,295 peptides   86.3% purity
	 * </pre>
	 *
	 * +15.3% PSMs and +5.4% peptides at the same FDR, with purity slightly
	 * better rather than worse. The extra ranks are not extra identifications -
	 * the spectrum collapse still keeps exactly one PSM per spectrum - they are a
	 * choice for the classifier, which overrules the engine's top hit for 15.0%
	 * of the final set. 9,195 peptides in the result are ones rank 1 never
	 * proposed, 87% of them on a true genome.
	 * <p>
	 * Round 1 stays at 1 ({@code round1OutputReportTopN}): it has no classifier,
	 * so nothing there can promote a lower-ranked candidate.
	 */
	// REVIEWED 2026-09-17 (Claude session, source-verified against the javadoc above): rank
	// policy re-confirmed, no open question. See Claude outputs/opensearch_optimization_review.md.
	private int round2OutputReportTopN = 5;
	private int maxIsotopeOffset = 3;
	/** Stage-1 folding is bounded separately; +2 is the last measured offset. */
	private int round1MaxIsotopeOffset = 2;
	private int minIsotopeFamilySatellites = 2;
	private double isotopeMaxAbundanceIncrease = 1.25;

	// --- stage 1: expected mass shifts (see ExpectedPeakSelector) -----------
	/** Off restores the previous zero-family-only stage-1 behaviour exactly. */
	private boolean expectedPeaksEnabled = true;
	/** Half-width of each shift's window (Da); 8 mDa separates +42.011 from +42.047. */
	private double expectedPeakToleranceDa = 0.008;
	/** Minimum target PSMs in a shift's window for it to be testable at all. */
	private int expectedPeakMinTargets = 20;
	/** Maximum raw decoys/targets in the window, before any cut. */
	private double expectedPeakMaxWindowFdr = 0.50;
	/** Minimum share of a shift's accepted peptides also seen unmodified this run. */
	private double expectedPeakMinCorroboration = 0.50;
	/** Cap on admitted targets, as a fraction of the zero family's accepted targets. */
	private double expectedPeakMaxAdmittedFraction = 0.25;

	// --- stage 1: rescoring (see Stage1Rescorer / MSBoosterHeadlessRunner) ---
	/**
	 * Replace stage 1's engine-score cut with a cross-validated semi-supervised
	 * rescoring. Off restores the previous behaviour exactly.
	 */
	private boolean stage1RescoringEnabled = true;
	/**
	 * Add MSBooster's predicted-spectrum and retention-time features to that
	 * rescoring. Requires a configured FragPipe (for the jar) and DIA-NN (as the
	 * prediction backend); when either is missing the rescorer runs on search
	 * features alone rather than failing.
	 */
	private boolean stage1MsboosterEnabled = true;
	/** Heap for the MSBooster JVM. Prediction is the memory-hungry part. */
	private String msboosterMaxHeap = "16g";
	/**
	 * Tool paths for the MSBooster step, set by the caller.
	 * <p>
	 * Resolved by the GUI/CLI/HPC entry point from the configured resources and
	 * placed here, exactly as {@code sageExecutable} is, so every entry point
	 * takes the executables from the same configuration and core never reaches
	 * into the Resources model itself.
	 */
	private String msboosterJar = "";
	private String diannExecutable = "";
	private String unimodObo = "";

	// --- stage 2: delta-mass evidence and substitution analysis ------------
	/**
	 * Reject fitted peaks whose members show no shifted-fragment support and fold
	 * their PSMs onto the zero family, rather than reporting them as
	 * modifications. See {@code DeltaMassEvidenceGate}. On by default: measured on
	 * 12mix (2026-09-14) this is what separates the +3.9..+4.3 and +4.8..+5.4 Da
	 * precursor artifacts - 94-96% of their PSMs gain nothing from placing the
	 * delta on fragment ions - from real modifications at 8-41%.
	 */
	private boolean deltaMassEvidenceGateEnabled = true;
	/** Share of a peak's members with no shifted-fragment gain that condemns it. */
	private double deltaMassMaxNoGainFraction = 0.80;
	/**
	 * Flag peaks whose delta is explained by a single amino-acid substitution and
	 * keep them out of the modification verdict. See {@code SubstitutionReassigner}.
	 * On by default: the +/-14.0157 Da peaks passed the per-peak local-FDR walk at
	 * 0.94%/0.96% decoy while their true error against the twelve spiked genomes
	 * was 18.35%/17.85%, because target-decoy cannot model a wrong peptide that is
	 * itself a target.
	 */
	private boolean substitutionAnalysisEnabled = true;
	/** Tolerance for matching a peak centroid to a substitution mass (Da). */
	private double substitutionMassToleranceDa = 0.005;
	/** Share of a peak's members yielding a localization-consistent candidate. */
	private double substitutionMinFraction = 0.50;
	/**
	 * Share of a peak's candidate sequences that must exist in the searched
	 * database before the peak is called a substitution rather than a
	 * modification. Without it, deamidation, +40.0313 and carbamidomethyl are all
	 * condemned - their colliding substitutions act on the same residues they do,
	 * so localization cannot separate them. See
	 * {@code SubstitutionReassigner.Config.minInDatabaseFraction}.
	 */
	private double substitutionMinInDatabaseFraction = 0.05;
	/**
	 * Call reproducible single-residue sequence variants from the substitution
	 * PSMs whose candidate is absent from every searched database, and write
	 * {@code sequence_variants.tsv} / {@code sequence_variants.fasta}. Reporting
	 * only: no PSM, q-value or quantification changes. See
	 * {@code SequenceVariantCaller}.
	 */
	private boolean sequenceVariantCallingEnabled = true;
	/**
	 * Distinct samples a variant site must be seen in to be called. On 12mix
	 * (9 files) sites seen in >= 3 files had a 3% two-nucleotide rate against
	 * 10-13% for sites seen in 1-2 files.
	 */
	private int sequenceVariantMinSamples = 3;
	/**
	 * Pangenome variant check (UHGP-100): report-only add-on run after the
	 * taxonomy step. Its MS2-undecidable test uses
	 * {@link #substitutionMassToleranceDa}.
	 */
	private final uk.ac.quadram.metapilot.core.variant.PangenomeVariantSettings pangenomeVariantSettings =
			new uk.ac.quadram.metapilot.core.variant.PangenomeVariantSettings();
	/**
	 * MS2 form ranking of the pangenome check's candidates: substitution
	 * against same-mass modification and other substitutions, within each
	 * spectrum, with DIA-NN predicted intensities ({@link #diannExecutable}).
	 */
	private final uk.ac.quadram.metapilot.core.variant.VariantFormRankingSettings variantFormRankingSettings =
			new uk.ac.quadram.metapilot.core.variant.VariantFormRankingSettings();

	// --- stage 2: Crystal-C precursor refinement ---------------------------
	/**
	 * Refine each round-two pepXML with Crystal-C before it reaches the pooled
	 * validation.
	 * <p>
	 * On by default because without it the workflow reports precursor-assignment
	 * artifacts as modifications. Measured on 12mix (2026-09-14): MetaPilot's
	 * open search and FragPipe's find the same number of rank-1 unmodified
	 * target PSMs in the raw pepXML - 352,484 against 352,360 - but FragPipe runs
	 * Crystal-C and MetaPilot did not, so its zero-delta count went to 591,641
	 * while MetaPilot's validated output left 316,120 of 785,180 target PSMs
	 * sitting at 0.05-5.5 Da and spent 335 of its 500 fitted peaks inside
	 * +/-5.5 Da. Off restores the previous behaviour exactly.
	 */
	private boolean crystalcEnabled = true;
	/**
	 * Heap for the Crystal-C JVM. It holds a whole run's MS1 scans while it
	 * re-picks monoisotopic peaks, so this is the memory-hungry part; FragPipe
	 * hands it the machine's whole budget and there is a standing upstream issue
	 * about it running out on large files.
	 */
	private String crystalcMaxHeap = "32g";
	/**
	 * Path to Crystal-C's jar, set by the caller from the configured FragPipe
	 * installation, exactly as {@code msboosterJar} is.
	 */
	private String crystalcJar = "";
	/**
	 * Crystal-C's own numeric settings, at FragPipe's bundled defaults. FragPipe
	 * overrides only the four paths and exposes none of these, so leaving them
	 * alone is what makes a MetaPilot run refine identically to the FragPipe run
	 * the behaviour above was measured against. Change them only with a
	 * measurement.
	 */
	private int crystalcMinCharge = 1;
	private int crystalcMaxCharge = 6;
	private int crystalcIsotopeNumber = 3;
	private double crystalcPrecursorMassPpm = 20;
	private double crystalcIsolationWindow = 0.7;

	// --- stage 2: validation -----------------------------------------------
	private int classifierFolds = 3;
	private int classifierMaxIterations = 10;
	private double classifierTrainQValue = 0.01;
	private double localFdr = 0.01;
	/** Minimum R-squared for a fitted delta-mass peak to be kept. */
	private double peakMinR2 = 0.9;

	// --- taxonomy ----------------------------------------------------------
	private boolean eukaryota = false;
	private boolean bacteria = false;
	private boolean archaea = false;
	private boolean viruses = false;
	private boolean viroids = false;
	private boolean unclassified = false;
	private boolean other = false;
	private String ignoreRank = "Species";
	private boolean ignoreEnvironmentalSamples = false;

	/** MetadataData format: [sampleName, labelType, age, gender, ...]. */
	private List<String[]> metadataTableData = new ArrayList<>();

	public DDAMagOpenSearchParameters() {
		super("DDA_MAG_OPENSEARCH", new MetaData(), "", "", false, "", "HCD-FTMS",
				Runtime.getRuntime().availableProcessors(), Runtime.getRuntime().availableProcessors(), true);
		applyDefaults();
	}

	public DDAMagOpenSearchParameters(String analysisType, MetaData metadata, String result, String microDb,
			boolean appendHostDb, String hostDB, String ms2ScanMode, int threadCount) {
		super(analysisType, metadata, result, microDb, appendHostDb, hostDB, ms2ScanMode, threadCount, threadCount,
				true);
		applyDefaults();
	}

	private void applyDefaults() {
		setEnzyme(MSFraggerHeadlessRunner.DEFAULT_ENZYME);
		setMissCleavages(2);
		setDigestMode("Specific");
		setMs2ScanMode(MetaConstants.HCD_FTMS);
		setThreadCount(Runtime.getRuntime().availableProcessors());
		setQuantitativeMode(MetaConstants.labelFree);
		setReduceDb(true);
		setSpectraCluster(false);
		setFixMods(new String[] { "Carbamidomethyl[C] 0" });
		// Open search finds shifted masses through the wide precursor window
		// itself, so an enumerated variable-mod list is not how modifications are
		// discovered here. Only the near-ubiquitous ones are declared.
		setVariMods(new String[] { "Oxidation[M] 0" });
	}

	/** Fitter configuration derived from these parameters. */
	public DeltaMassPeakFitter.Config toFitterConfig() {
		DeltaMassPeakFitter.Config cfg = new DeltaMassPeakFitter.Config();
		cfg.windowLower = precursorMassLower;
		cfg.windowUpper = precursorMassUpper;
		cfg.membershipSigmas = zeroPeakSigmas;
		cfg.membershipFloorDa = membershipFloorDa;
		cfg.zeroWindow = zeroPeakSearchWindow;
		cfg.minR2 = peakMinR2;
		cfg.minPeakArea = peakMinArea;
		cfg.maxPeaksToFit = maxPeaksToFit;
		cfg.maxFitIterations = maxFitIterations;
		cfg.minFitBins = peakMinFitBins;
		cfg.minSigmaBins = peakMinSigmaBins;
		cfg.minPeakSnr = peakMinSnr;
		cfg.minPeakSeparationBins = peakMinSeparationBins;
		cfg.minPeakSeparationSigmas = peakMinSeparationSigmas;
		cfg.useKde = useKdePeakDiscovery;
		cfg.kdeBandwidthDa = kdeBandwidthDa;
		cfg.kdeScaleMultipliers = kdeScaleMultipliers;
		cfg.kdeMinPersistentScales = kdeMinPersistentScales;
		cfg.kdeMaxMembershipDa = kdeMaxMembershipDa;
		cfg.kdeMinReproducibleFiles = kdeMinReproducibleFiles;
		cfg.kdeMinFileSupport = kdeMinFileSupport;
		cfg.kdeMaxModes = kdeMaxModes;
		cfg.kdeMaxSigmaMultiple = kdeMaxSigmaMultiple;
		cfg.kdeMinProminenceRatio = kdeMinProminenceRatio;
		cfg.maxIsotopeOffset = maxIsotopeOffset;
		cfg.round1MaxIsotopeOffset = round1MaxIsotopeOffset;
		cfg.minIsotopeFamilySatellites = minIsotopeFamilySatellites;
		cfg.isotopeMaxAbundanceIncrease = isotopeMaxAbundanceIncrease;
		return cfg;
	}

	public double getMembershipFloorDa() {
		return membershipFloorDa;
	}

	public void setMembershipFloorDa(double membershipFloorDa) {
		this.membershipFloorDa = membershipFloorDa;
	}

	public int getPeakMinArea() {
		return peakMinArea;
	}

	public void setPeakMinArea(int peakMinArea) {
		this.peakMinArea = peakMinArea;
	}

	public int getMaxPeaksToFit() {
		return maxPeaksToFit;
	}

	public void setMaxPeaksToFit(int maxPeaksToFit) {
		this.maxPeaksToFit = maxPeaksToFit;
	}

	public int getMaxFitIterations() {
		return maxFitIterations;
	}

	public void setMaxFitIterations(int maxFitIterations) {
		this.maxFitIterations = maxFitIterations;
	}

	public int getPeakMinFitBins() {
		return peakMinFitBins;
	}

	public void setPeakMinFitBins(int peakMinFitBins) {
		this.peakMinFitBins = peakMinFitBins;
	}

	public double getPeakMinSigmaBins() {
		return peakMinSigmaBins;
	}

	public void setPeakMinSigmaBins(double peakMinSigmaBins) {
		this.peakMinSigmaBins = peakMinSigmaBins;
	}

	public double getPeakMinSnr() {
		return peakMinSnr;
	}

	public void setPeakMinSnr(double peakMinSnr) {
		this.peakMinSnr = peakMinSnr;
	}

	public double getPeakMinSeparationBins() {
		return peakMinSeparationBins;
	}

	public void setPeakMinSeparationBins(double peakMinSeparationBins) {
		this.peakMinSeparationBins = peakMinSeparationBins;
	}

	public double getPeakMinSeparationSigmas() {
		return peakMinSeparationSigmas;
	}

	public void setPeakMinSeparationSigmas(double peakMinSeparationSigmas) {
		this.peakMinSeparationSigmas = peakMinSeparationSigmas;
	}

	public boolean isUseKdePeakDiscovery() {
		return useKdePeakDiscovery;
	}

	public void setUseKdePeakDiscovery(boolean useKdePeakDiscovery) {
		this.useKdePeakDiscovery = useKdePeakDiscovery;
	}

	public double getKdeBandwidthDa() {
		return kdeBandwidthDa;
	}

	public void setKdeBandwidthDa(double kdeBandwidthDa) {
		this.kdeBandwidthDa = kdeBandwidthDa;
	}

	public String getKdeScaleMultipliers() {
		return kdeScaleMultipliers;
	}

	public void setKdeScaleMultipliers(String kdeScaleMultipliers) {
		this.kdeScaleMultipliers = kdeScaleMultipliers;
	}

	public int getKdeMinPersistentScales() {
		return kdeMinPersistentScales;
	}

	public void setKdeMinPersistentScales(int kdeMinPersistentScales) {
		this.kdeMinPersistentScales = kdeMinPersistentScales;
	}

	public double getKdeMaxMembershipDa() {
		return kdeMaxMembershipDa;
	}

	public void setKdeMaxMembershipDa(double kdeMaxMembershipDa) {
		this.kdeMaxMembershipDa = kdeMaxMembershipDa;
	}

	public int getKdeMinReproducibleFiles() {
		return kdeMinReproducibleFiles;
	}

	public void setKdeMinReproducibleFiles(int kdeMinReproducibleFiles) {
		this.kdeMinReproducibleFiles = kdeMinReproducibleFiles;
	}

	public double getKdeMinFileSupport() {
		return kdeMinFileSupport;
	}

	public void setKdeMinFileSupport(double kdeMinFileSupport) {
		this.kdeMinFileSupport = kdeMinFileSupport;
	}

	public int getKdeMaxModes() {
		return kdeMaxModes;
	}

	public double getKdeMinProminenceRatio() {
		return kdeMinProminenceRatio;
	}

	public void setKdeMinProminenceRatio(double kdeMinProminenceRatio) {
		this.kdeMinProminenceRatio = kdeMinProminenceRatio;
	}

	public double getKdeMaxSigmaMultiple() {
		return kdeMaxSigmaMultiple;
	}

	public void setKdeMaxSigmaMultiple(double kdeMaxSigmaMultiple) {
		this.kdeMaxSigmaMultiple = kdeMaxSigmaMultiple;
	}

	public boolean isCombinedPsmInstances() {
		return combinedPsmInstances;
	}

	public void setCombinedPsmInstances(boolean combinedPsmInstances) {
		this.combinedPsmInstances = combinedPsmInstances;
	}

	/** True when this run should merge engine rows into combined candidates before rescoring. */
	public boolean usesCombinedPsmInstances() {
		return combinedPsmInstances && usesSage();
	}

	public boolean isRestrictPositivesToValidatedPeaks() {
		return restrictPositivesToValidatedPeaks;
	}

	public void setRestrictPositivesToValidatedPeaks(boolean restrictPositivesToValidatedPeaks) {
		this.restrictPositivesToValidatedPeaks = restrictPositivesToValidatedPeaks;
	}

	public int getRound2OutputReportTopN() {
		return round2OutputReportTopN;
	}

	public void setRound2OutputReportTopN(int round2OutputReportTopN) {
		this.round2OutputReportTopN = round2OutputReportTopN;
	}

	public void setKdeMaxModes(int kdeMaxModes) {
		this.kdeMaxModes = kdeMaxModes;
	}

	public int getMaxIsotopeOffset() {
		return maxIsotopeOffset;
	}

	public void setMaxIsotopeOffset(int maxIsotopeOffset) {
		this.maxIsotopeOffset = maxIsotopeOffset;
	}

	public int getRound1MaxIsotopeOffset() {
		return round1MaxIsotopeOffset;
	}

	public void setRound1MaxIsotopeOffset(int round1MaxIsotopeOffset) {
		this.round1MaxIsotopeOffset = round1MaxIsotopeOffset;
	}

	public int getMinIsotopeFamilySatellites() {
		return minIsotopeFamilySatellites;
	}

	public void setMinIsotopeFamilySatellites(int minIsotopeFamilySatellites) {
		this.minIsotopeFamilySatellites = minIsotopeFamilySatellites;
	}

	public double getIsotopeMaxAbundanceIncrease() {
		return isotopeMaxAbundanceIncrease;
	}

	public void setIsotopeMaxAbundanceIncrease(double isotopeMaxAbundanceIncrease) {
		this.isotopeMaxAbundanceIncrease = isotopeMaxAbundanceIncrease;
	}

	/** Validator configuration derived from these parameters. */
	public ModificationValidator.Config toValidatorConfig() {
		ModificationValidator.Config cfg = new ModificationValidator.Config();
		cfg.localFdr = localFdr;
		cfg.unexplainedFdr = localFdr;
		return cfg;
	}

	// --- accessors ---------------------------------------------------------

	public List<String[]> getExperimentData() {
		return experimentData;
	}

	public void setExperimentData(List<String[]> experimentData) {
		this.experimentData = experimentData == null ? new ArrayList<>() : experimentData;
	}

	public MagDbItem getMagDbItem() {
		return magDbItem;
	}

	public void setMagDbItem(MagDbItem magDbItem) {
		this.magDbItem = magDbItem;
	}

	public String getMagDbCatalogueId() {
		return magDbCatalogueId;
	}

	public void setMagDbCatalogueId(String magDbCatalogueId) {
		this.magDbCatalogueId = magDbCatalogueId;
	}

	public String getMagDbIdentifier() {
		return magDbIdentifier;
	}

	public void setMagDbIdentifier(String magDbIdentifier) {
		this.magDbIdentifier = magDbIdentifier;
	}

	public String getMagDbVersion() {
		return magDbVersion;
	}

	public void setMagDbVersion(String magDbVersion) {
		this.magDbVersion = magDbVersion;
	}

	public String getGenomeSelectorMode() {
		return genomeSelectorMode;
	}

	public void setGenomeSelectorMode(String genomeSelectorMode) {
		this.genomeSelectorMode = genomeSelectorMode;
	}

	public String getGenomeSelectionStrategy() {
		return genomeSelectionStrategy;
	}

	public void setGenomeSelectionStrategy(String genomeSelectionStrategy) {
		this.genomeSelectionStrategy = genomeSelectionStrategy;
	}

	public String getSearchEngines() {
		return searchEngines;
	}

	public void setSearchEngines(String searchEngines) {
		this.searchEngines = searchEngines == null || searchEngines.isBlank() ? ENGINES_MSFRAGGER : searchEngines;
	}

	public boolean isMsfraggerOnly() {
		return ENGINES_MSFRAGGER.equals(searchEngines);
	}

	/** True when the selected mode includes Sage as a second engine. */
	public boolean usesSage() {
		return ENGINES_MSFRAGGER_AND_SAGE.equals(searchEngines);
	}

	/**
	 * Path to the Sage executable, supplied by the caller from configured
	 * resources. Kept on the parameters so GUI, CLI and HPC entry points resolve
	 * it the same way, as they do for MSFragger.
	 */
	private String sageExecutable = "";

	public String getSageExecutable() {
		return sageExecutable;
	}

	public void setSageExecutable(String sageExecutable) {
		this.sageExecutable = sageExecutable == null ? "" : sageExecutable;
	}

	/**
	 * Path to ProteoWizard's msconvert, supplied by the caller from configured
	 * resources like {@link #getSageExecutable()}. Needed only when Thermo
	 * {@code .raw} files are loaded: they are converted to mzML before the
	 * search.
	 */
	private String msconvertExecutable = "";

	public String getMsconvertExecutable() {
		return msconvertExecutable;
	}

	public void setMsconvertExecutable(String msconvertExecutable) {
		this.msconvertExecutable = msconvertExecutable == null ? "" : msconvertExecutable;
	}

	/**
	 * Whether an engine mode can actually run today.
	 * <p>
	 * Single source of truth for what is implemented. The UI, the pre-flight
	 * check and the headless runner all consult this, so enabling Sage means
	 * changing one method rather than finding every place that happened to
	 * hardcode the same condition.
	 */
	public static boolean isEngineModeAvailable(String engineMode) {
		return ENGINES_MSFRAGGER.equals(engineMode) || ENGINES_MSFRAGGER_AND_SAGE.equals(engineMode);
	}

	/** Why {@code engineMode} cannot run, or null when it can. */
	public static String engineModeUnavailableReason(String engineMode) {
		if (isEngineModeAvailable(engineMode)) {
			return null;
		}
		return "Engine mode '" + engineMode + "' is not recognised. Available modes: "
				+ ENGINES_MSFRAGGER + ", " + ENGINES_MSFRAGGER_AND_SAGE + ".";
	}

	/**
	 * Coerces an unavailable mode back to the only runnable one.
	 * <p>
	 * Applied when loading saved parameters: a file written while the mode was
	 * selectable, or hand-edited, must not put the workflow into a state that
	 * only fails at run time.
	 */
	public static String coerceEngineMode(String engineMode) {
		return isEngineModeAvailable(engineMode) ? engineMode : ENGINES_MSFRAGGER;
	}

	public int getParallelTasks() {
		return parallelTasks;
	}

	public void setParallelTasks(int parallelTasks) {
		this.parallelTasks = parallelTasks;
	}

	public String getMsfraggerMaxHeap() {
		return msfraggerMaxHeap;
	}

	public void setMsfraggerMaxHeap(String msfraggerMaxHeap) {
		this.msfraggerMaxHeap = msfraggerMaxHeap;
	}

	public double getPrecursorMassLower() {
		return precursorMassLower;
	}

	public void setPrecursorMassLower(double precursorMassLower) {
		this.precursorMassLower = precursorMassLower;
	}

	public double getPrecursorMassUpper() {
		return precursorMassUpper;
	}

	public void setPrecursorMassUpper(double precursorMassUpper) {
		this.precursorMassUpper = precursorMassUpper;
	}

	public int getOutputReportTopN() {
		return outputReportTopN;
	}

	public void setOutputReportTopN(int outputReportTopN) {
		this.outputReportTopN = outputReportTopN;
	}

	public double getZeroPeakSigmas() {
		return zeroPeakSigmas;
	}

	public void setZeroPeakSigmas(double zeroPeakSigmas) {
		this.zeroPeakSigmas = zeroPeakSigmas;
	}

	public double getZeroPeakSearchWindow() {
		return zeroPeakSearchWindow;
	}

	public void setZeroPeakSearchWindow(double zeroPeakSearchWindow) {
		this.zeroPeakSearchWindow = zeroPeakSearchWindow;
	}

	public int getClassifierFolds() {
		return classifierFolds;
	}

	public void setClassifierFolds(int classifierFolds) {
		this.classifierFolds = classifierFolds;
	}

	public int getClassifierMaxIterations() {
		return classifierMaxIterations;
	}

	public void setClassifierMaxIterations(int classifierMaxIterations) {
		this.classifierMaxIterations = classifierMaxIterations;
	}

	public double getClassifierTrainQValue() {
		return classifierTrainQValue;
	}

	public void setClassifierTrainQValue(double classifierTrainQValue) {
		this.classifierTrainQValue = classifierTrainQValue;
	}







	public double getLocalFdr() {
		return localFdr;
	}

	public void setLocalFdr(double localFdr) {
		this.localFdr = localFdr;
	}

	public double getPeakMinR2() {
		return peakMinR2;
	}

	public void setPeakMinR2(double peakMinR2) {
		this.peakMinR2 = peakMinR2;
	}

	public boolean isEukaryota() {
		return eukaryota;
	}

	public void setEukaryota(boolean eukaryota) {
		this.eukaryota = eukaryota;
	}

	public boolean isBacteria() {
		return bacteria;
	}

	public void setBacteria(boolean bacteria) {
		this.bacteria = bacteria;
	}

	public boolean isArchaea() {
		return archaea;
	}

	public void setArchaea(boolean archaea) {
		this.archaea = archaea;
	}

	public boolean isViruses() {
		return viruses;
	}

	public void setViruses(boolean viruses) {
		this.viruses = viruses;
	}

	public boolean isViroids() {
		return viroids;
	}

	public void setViroids(boolean viroids) {
		this.viroids = viroids;
	}

	public boolean isUnclassified() {
		return unclassified;
	}

	public void setUnclassified(boolean unclassified) {
		this.unclassified = unclassified;
	}

	public boolean isOther() {
		return other;
	}

	public void setOther(boolean other) {
		this.other = other;
	}

	public String getIgnoreRank() {
		return ignoreRank;
	}

	public void setIgnoreRank(String ignoreRank) {
		this.ignoreRank = ignoreRank;
	}

	public boolean isIgnoreEnvironmentalSamples() {
		return ignoreEnvironmentalSamples;
	}

	public void setIgnoreEnvironmentalSamples(boolean ignoreEnvironmentalSamples) {
		this.ignoreEnvironmentalSamples = ignoreEnvironmentalSamples;
	}

	public List<String[]> getMetadataTableData() {
		return metadataTableData;
	}

	public void setMetadataTableData(List<String[]> metadataTableData) {
		this.metadataTableData = metadataTableData == null ? new ArrayList<>() : metadataTableData;
	}

	public int getRound1OutputReportTopN() {
		return round1OutputReportTopN;
	}

	public void setRound1OutputReportTopN(int round1OutputReportTopN) {
		this.round1OutputReportTopN = round1OutputReportTopN;
	}

	public String getRound1SearchMode() {
		return round1SearchMode;
	}

	/** Anything other than {@code open} (case-insensitive) is taken as {@code closed}. */
	public void setRound1SearchMode(String round1SearchMode) {
		this.round1SearchMode = ROUND1_SEARCH_OPEN.equalsIgnoreCase(
				round1SearchMode == null ? "" : round1SearchMode.trim()) ? ROUND1_SEARCH_OPEN : ROUND1_SEARCH_CLOSED;
	}

	public boolean isRound1ClosedSearch() {
		return ROUND1_SEARCH_CLOSED.equals(round1SearchMode);
	}

	public double getRound1ClosedPrecursorPpm() {
		return round1ClosedPrecursorPpm;
	}

	public void setRound1ClosedPrecursorPpm(double round1ClosedPrecursorPpm) {
		this.round1ClosedPrecursorPpm = round1ClosedPrecursorPpm;
	}

	public boolean isClosedStageEnabled() {
		return closedStageEnabled;
	}

	public void setClosedStageEnabled(boolean closedStageEnabled) {
		this.closedStageEnabled = closedStageEnabled;
	}

	public double getClosedStagePrecursorPpm() {
		return closedStagePrecursorPpm;
	}

	public void setClosedStagePrecursorPpm(double closedStagePrecursorPpm) {
		this.closedStagePrecursorPpm = closedStagePrecursorPpm;
	}

	public boolean isExpectedPeaksEnabled() {
		return expectedPeaksEnabled;
	}

	public void setExpectedPeaksEnabled(boolean expectedPeaksEnabled) {
		this.expectedPeaksEnabled = expectedPeaksEnabled;
	}

	public double getExpectedPeakToleranceDa() {
		return expectedPeakToleranceDa;
	}

	public void setExpectedPeakToleranceDa(double expectedPeakToleranceDa) {
		this.expectedPeakToleranceDa = expectedPeakToleranceDa;
	}

	public int getExpectedPeakMinTargets() {
		return expectedPeakMinTargets;
	}

	public void setExpectedPeakMinTargets(int expectedPeakMinTargets) {
		this.expectedPeakMinTargets = expectedPeakMinTargets;
	}

	public double getExpectedPeakMaxWindowFdr() {
		return expectedPeakMaxWindowFdr;
	}

	public void setExpectedPeakMaxWindowFdr(double expectedPeakMaxWindowFdr) {
		this.expectedPeakMaxWindowFdr = expectedPeakMaxWindowFdr;
	}

	public double getExpectedPeakMinCorroboration() {
		return expectedPeakMinCorroboration;
	}

	public void setExpectedPeakMinCorroboration(double expectedPeakMinCorroboration) {
		this.expectedPeakMinCorroboration = expectedPeakMinCorroboration;
	}

	public double getExpectedPeakMaxAdmittedFraction() {
		return expectedPeakMaxAdmittedFraction;
	}

	public void setExpectedPeakMaxAdmittedFraction(double expectedPeakMaxAdmittedFraction) {
		this.expectedPeakMaxAdmittedFraction = expectedPeakMaxAdmittedFraction;
	}

	/**
	 * Stage-1 expected-shift criteria, in the same shape
	 * {@link #toFitterConfig()} builds the fitter's.
	 */
	public ExpectedPeakSelector.Config toExpectedPeakConfig() {
		ExpectedPeakSelector.Config cfg = new ExpectedPeakSelector.Config();
		cfg.enabled = expectedPeaksEnabled;
		cfg.toleranceDa = expectedPeakToleranceDa;
		cfg.minTargetsInWindow = expectedPeakMinTargets;
		cfg.maxWindowFdr = expectedPeakMaxWindowFdr;
		cfg.minCorroboration = expectedPeakMinCorroboration;
		cfg.maxAdmittedFraction = expectedPeakMaxAdmittedFraction;
		return cfg;
	}

	public boolean isStage1RescoringEnabled() {
		return stage1RescoringEnabled;
	}

	public void setStage1RescoringEnabled(boolean stage1RescoringEnabled) {
		this.stage1RescoringEnabled = stage1RescoringEnabled;
	}

	public boolean isStage1MsboosterEnabled() {
		return stage1MsboosterEnabled;
	}

	public void setStage1MsboosterEnabled(boolean stage1MsboosterEnabled) {
		this.stage1MsboosterEnabled = stage1MsboosterEnabled;
	}

	public String getMsboosterMaxHeap() {
		return msboosterMaxHeap;
	}

	public void setMsboosterMaxHeap(String msboosterMaxHeap) {
		this.msboosterMaxHeap = msboosterMaxHeap;
	}

	/** Stage-1 rescoring configuration, in the shape {@link Stage1Rescorer} takes. */
	public Stage1Rescorer.Config toStage1RescorerConfig() {
		Stage1Rescorer.Config cfg = new Stage1Rescorer.Config();
		cfg.enabled = stage1RescoringEnabled;
		cfg.fdrThreshold = UnmodifiedPsmSelector.DEFAULT_FDR;
		cfg.classifier.folds = classifierFolds;
		cfg.classifier.maxIterations = classifierMaxIterations;
		cfg.classifier.trainQValue = classifierTrainQValue;
		return cfg;
	}

	public String getMsboosterJarPath() {
		return msboosterJar;
	}

	public void setMsboosterJarPath(String msboosterJar) {
		this.msboosterJar = msboosterJar == null ? "" : msboosterJar;
	}

	/** The configured MSBooster jar, or null when it is unset or missing. */
	public java.io.File getMsboosterJar() {
		return asExistingFile(msboosterJar);
	}

	public String getDiannExecutablePath() {
		return diannExecutable;
	}

	public void setDiannExecutablePath(String diannExecutable) {
		this.diannExecutable = diannExecutable == null ? "" : diannExecutable;
	}

	/** The configured DIA-NN executable, or null when it is unset or missing. */
	public java.io.File getDiannExecutable() {
		return asExistingFile(diannExecutable);
	}

	public String getUnimodOboPath() {
		return unimodObo;
	}

	public void setUnimodOboPath(String unimodObo) {
		this.unimodObo = unimodObo == null ? "" : unimodObo;
	}

	/** {@code unimod.obo}, or null; MSBooster falls back to its own defaults. */
	public java.io.File getUnimodObo() {
		return asExistingFile(unimodObo);
	}

	public boolean isCrystalcEnabled() {
		return crystalcEnabled;
	}

	public void setCrystalcEnabled(boolean crystalcEnabled) {
		this.crystalcEnabled = crystalcEnabled;
	}

	public String getCrystalcMaxHeap() {
		return crystalcMaxHeap;
	}

	public void setCrystalcMaxHeap(String crystalcMaxHeap) {
		this.crystalcMaxHeap = crystalcMaxHeap;
	}

	public String getCrystalcJarPath() {
		return crystalcJar;
	}

	public void setCrystalcJarPath(String crystalcJar) {
		this.crystalcJar = crystalcJar == null ? "" : crystalcJar;
	}

	/** The configured Crystal-C jar, or null when it is unset or missing. */
	public java.io.File getCrystalcJar() {
		return asExistingFile(crystalcJar);
	}

	public int getCrystalcMinCharge() {
		return crystalcMinCharge;
	}

	public void setCrystalcMinCharge(int crystalcMinCharge) {
		this.crystalcMinCharge = crystalcMinCharge;
	}

	public int getCrystalcMaxCharge() {
		return crystalcMaxCharge;
	}

	public void setCrystalcMaxCharge(int crystalcMaxCharge) {
		this.crystalcMaxCharge = crystalcMaxCharge;
	}

	public int getCrystalcIsotopeNumber() {
		return crystalcIsotopeNumber;
	}

	public void setCrystalcIsotopeNumber(int crystalcIsotopeNumber) {
		this.crystalcIsotopeNumber = crystalcIsotopeNumber;
	}

	public double getCrystalcPrecursorMassPpm() {
		return crystalcPrecursorMassPpm;
	}

	public void setCrystalcPrecursorMassPpm(double crystalcPrecursorMassPpm) {
		this.crystalcPrecursorMassPpm = crystalcPrecursorMassPpm;
	}

	public double getCrystalcIsolationWindow() {
		return crystalcIsolationWindow;
	}

	public void setCrystalcIsolationWindow(double crystalcIsolationWindow) {
		this.crystalcIsolationWindow = crystalcIsolationWindow;
	}

	public boolean isDeltaMassEvidenceGateEnabled() {
		return deltaMassEvidenceGateEnabled;
	}

	public void setDeltaMassEvidenceGateEnabled(boolean deltaMassEvidenceGateEnabled) {
		this.deltaMassEvidenceGateEnabled = deltaMassEvidenceGateEnabled;
	}

	public double getDeltaMassMaxNoGainFraction() {
		return deltaMassMaxNoGainFraction;
	}

	public void setDeltaMassMaxNoGainFraction(double deltaMassMaxNoGainFraction) {
		this.deltaMassMaxNoGainFraction = deltaMassMaxNoGainFraction;
	}

	public boolean isSubstitutionAnalysisEnabled() {
		return substitutionAnalysisEnabled;
	}

	public void setSubstitutionAnalysisEnabled(boolean substitutionAnalysisEnabled) {
		this.substitutionAnalysisEnabled = substitutionAnalysisEnabled;
	}

	public double getSubstitutionMassToleranceDa() {
		return substitutionMassToleranceDa;
	}

	public void setSubstitutionMassToleranceDa(double substitutionMassToleranceDa) {
		this.substitutionMassToleranceDa = substitutionMassToleranceDa;
	}

	public double getSubstitutionMinFraction() {
		return substitutionMinFraction;
	}

	public void setSubstitutionMinFraction(double substitutionMinFraction) {
		this.substitutionMinFraction = substitutionMinFraction;
	}

	public double getSubstitutionMinInDatabaseFraction() {
		return substitutionMinInDatabaseFraction;
	}

	public void setSubstitutionMinInDatabaseFraction(double substitutionMinInDatabaseFraction) {
		this.substitutionMinInDatabaseFraction = substitutionMinInDatabaseFraction;
	}

	/** Config for {@code DeltaMassEvidenceGate}, in the shape it takes. */
	public uk.ac.quadram.metapilot.core.opensearch.DeltaMassEvidenceGate.Config toDeltaMassEvidenceConfig() {
		uk.ac.quadram.metapilot.core.opensearch.DeltaMassEvidenceGate.Config cfg =
				new uk.ac.quadram.metapilot.core.opensearch.DeltaMassEvidenceGate.Config();
		cfg.enabled = deltaMassEvidenceGateEnabled;
		cfg.maxNoGainFraction = deltaMassMaxNoGainFraction;
		cfg.zeroWindowDa = toFitterConfig().zeroWindow;
		return cfg;
	}

	/** Config for {@code SubstitutionReassigner}, in the shape it takes. */
	public uk.ac.quadram.metapilot.core.opensearch.SubstitutionReassigner.Config toSubstitutionConfig() {
		uk.ac.quadram.metapilot.core.opensearch.SubstitutionReassigner.Config cfg =
				new uk.ac.quadram.metapilot.core.opensearch.SubstitutionReassigner.Config();
		cfg.enabled = substitutionAnalysisEnabled;
		cfg.massToleranceDa = substitutionMassToleranceDa;
		cfg.minSubstitutionFraction = substitutionMinFraction;
		cfg.minInDatabaseFraction = substitutionMinInDatabaseFraction;
		return cfg;
	}

	public boolean isSequenceVariantCallingEnabled() {
		return sequenceVariantCallingEnabled;
	}

	public void setSequenceVariantCallingEnabled(boolean sequenceVariantCallingEnabled) {
		this.sequenceVariantCallingEnabled = sequenceVariantCallingEnabled;
	}

	public int getSequenceVariantMinSamples() {
		return sequenceVariantMinSamples;
	}

	public void setSequenceVariantMinSamples(int sequenceVariantMinSamples) {
		this.sequenceVariantMinSamples = sequenceVariantMinSamples;
	}

	/** Config for {@code SequenceVariantCaller}, in the shape it takes. */
	public uk.ac.quadram.metapilot.core.opensearch.SequenceVariantCaller.Config toSequenceVariantConfig() {
		uk.ac.quadram.metapilot.core.opensearch.SequenceVariantCaller.Config cfg =
				new uk.ac.quadram.metapilot.core.opensearch.SequenceVariantCaller.Config();
		// Needs the substitution verdicts; without them there is nothing to call.
		cfg.enabled = sequenceVariantCallingEnabled && substitutionAnalysisEnabled;
		cfg.minSamples = sequenceVariantMinSamples;
		cfg.isobaricToleranceDa = substitutionMassToleranceDa;
		return cfg;
	}

	public uk.ac.quadram.metapilot.core.variant.PangenomeVariantSettings getPangenomeVariantSettings() {
		return pangenomeVariantSettings;
	}

	public uk.ac.quadram.metapilot.core.variant.VariantFormRankingSettings getVariantFormRankingSettings() {
		return variantFormRankingSettings;
	}

	private static java.io.File asExistingFile(String path) {
		if (path == null || path.isBlank()) {
			return null;
		}
		java.io.File f = new java.io.File(path.trim());
		return f.exists() ? f : null;
	}
}
