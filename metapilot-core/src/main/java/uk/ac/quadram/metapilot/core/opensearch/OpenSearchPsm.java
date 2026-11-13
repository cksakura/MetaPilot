package uk.ac.quadram.metapilot.core.opensearch;

import java.util.Collections;
import java.util.LinkedHashMap;
import java.util.Locale;
import java.util.Map;

/**
 * One PSM candidate from a standalone MSFragger open search.
 * <p>
 * Holds the raw, unfiltered identification exactly as MSFragger reported it —
 * including candidates beyond rank 1 — plus the delta-mass peak context added
 * by {@link DeltaMassPeakFitter} and the score assigned by the validation
 * classifier. Nothing here is filtered: every field is either evidence for the
 * classifier or output from it.
 * <p>
 * Deliberately distinct from {@code dbSearch.msFragger.MSFraggerPsm}, which
 * carries the b/y ion-intensity fields the old MetaLab implementation derived
 * by re-matching raw MGF spectra by hand. Those are superseded here by
 * MSFragger's own localization-aware scores.
 */
public class OpenSearchPsm {

	// --- spectrum identity -------------------------------------------------
	private String fileName;
	private int scan;
	private int charge;
	private double precursorNeutralMass;
	private double retentionTime;

	/** 1 = best-scoring candidate for this spectrum; > 1 are the runners-up. */
	private int hitRank;

	// --- peptide -----------------------------------------------------------
	private String peptide;
	private String prevAA;
	private String nextAA;
	private double calcNeutralPepMass;
	private double massDiff;
	private int numTolTerm;
	private int numMissedCleavages;
	private int numMatchedIons;
	private int totNumIons;

	/**
	 * Which search engine produced this candidate.
	 * <p>
	 * Load-bearing, not bookkeeping: Sage cannot supply MSFragger's delta-mass
	 * localization block, so the two engines are scored by separate classifiers
	 * over their own feature sets and only pooled afterwards. This field is what
	 * partitions them. It is deliberately <em>not</em> a classifier feature — a
	 * model given the engine flag alongside a systematically-missing feature
	 * block would learn engine identity as a proxy for quality.
	 */
	private SearchEngine engine = SearchEngine.MSFRAGGER;

	public enum SearchEngine {
		MSFRAGGER, SAGE
	}

	public SearchEngine getEngine() {
		return engine;
	}

	public void setEngine(SearchEngine engine) {
		this.engine = engine == null ? SearchEngine.MSFRAGGER : engine;
	}

	// --- closed stage ------------------------------------------------------
	/**
	 * From the closed search of the closed-then-open identification. Such a PSM has no
	 * delta-mass peak ({@code peakId = -1}) and is not unexplained: its modifications are the
	 * search's variable modifications, in {@link #closedModifications}.
	 */
	private boolean closedStage;
	/**
	 * Closed stage only: variable modification mass per position, index 0 the peptide
	 * N-terminus and {@code i} residue {@code i} (1-based); null when there is none.
	 */
	private double[] closedModifications;

	public boolean isClosedStage() {
		return closedStage;
	}

	public void setClosedStage(boolean closedStage) {
		this.closedStage = closedStage;
	}

	public double[] getClosedModifications() {
		return closedModifications;
	}

	public void setClosedModifications(double[] closedModifications) {
		this.closedModifications = closedModifications;
	}

	public double getClosedModificationMass() {
		double total = 0.0;
		if (closedModifications != null) {
			for (double m : closedModifications) {
				total += m;
			}
		}
		return total;
	}

	public int getClosedModificationCount() {
		int n = 0;
		if (closedModifications != null) {
			for (double m : closedModifications) {
				if (m != 0.0) {
					n++;
				}
			}
		}
		return n;
	}

	/**
	 * The peptide with its variable modifications, {@code [+42.0106]-PEPM[+15.9949]K}: the
	 * N-terminal one before the sequence, a residue's after the residue. The bare peptide when
	 * unmodified, so a closed and an open unmodified PSM of one peptide share a Full Sequence.
	 */
	public String getClosedModifiedSequence() {
		if (peptide == null || getClosedModificationCount() == 0) {
			return peptide;
		}
		StringBuilder sb = new StringBuilder();
		if (closedModifications[0] != 0.0) {
			sb.append(String.format(Locale.ROOT, "[%+.4f]-", closedModifications[0]));
		}
		for (int i = 0; i < peptide.length(); i++) {
			sb.append(peptide.charAt(i));
			if (i + 1 < closedModifications.length && closedModifications[i + 1] != 0.0) {
				sb.append(String.format(Locale.ROOT, "[%+.4f]", closedModifications[i + 1]));
			}
		}
		return sb.toString();
	}

	// --- protein -----------------------------------------------------------
	/** pepXML's delegate protein, i.e. the {@code protein} attribute of search_hit. */
	private String protein;
	/**
	 * Every protein this peptide maps to: the delegate plus all
	 * {@code <alternative_protein>} entries.
	 * <p>
	 * Load-bearing for MAG genome selection, which counts how many genomes each
	 * peptide maps to. Keeping only the delegate would credit a marker peptide
	 * shared across N genomes to just one of them. The FragPipe path gets this
	 * for free because Philosopher's psm.tsv already flattens the delegate and
	 * "Mapped Proteins" into one list; bypassing Philosopher means doing that
	 * flattening here instead.
	 */
	private String[] allProteins = new String[0];
	private boolean target;

	// --- core scores -------------------------------------------------------
	private double hyperscore = Double.NaN;
	private double nextscore = Double.NaN;
	private double expect = Double.NaN;

	// --- delta-mass localization (MSFragger, localize_delta_mass = 1) -------
	/** Peptide sequence with probable delta-mass sites lowercased. */
	private String localizedSequence;
	private String deltaMassReassignment;
	private String substitutionVerdict;
	private String substitutionCandidate;
	private double scoreWithoutDeltaMass = Double.NaN;
	private double scoreWithDeltaMass = Double.NaN;
	private double secondBestScoreWithDeltaMass = Double.NaN;
	private double ionsWithoutDeltaMass = Double.NaN;
	private double ionsWithDeltaMass = Double.NaN;
	/** Every attribute from MSFragger's {@code ptm_result} element. */
	private Map<String, String> ptmResultAttributes = Collections.emptyMap();

	/**
	 * Every {@code search_score} element as reported, by name. Retained whole
	 * because MSFragger's pepXML score vocabulary is undocumented and version
	 * dependent — unrecognised entries stay available rather than being lost.
	 */
	private Map<String, String> searchScores = Collections.emptyMap();

	// --- delta-mass peak context (from DeltaMassPeakFitter) -----------------
	/** Raw KDE/Gaussian mode containing the observed precursor mass difference. */
	private int observedPeakId = -1;
	/** Canonical isotope-family root. This remains {@code -1} for unexplained PSMs. */
	private int peakId = -1;

	/**
	 * True when this PSM's peak passed the held-out gate. Distinct from
	 * {@link #isInPeak()}: a fitted peak is a location the discovery step
	 * proposed, a validated peak is one the held-out files then confirmed, and
	 * only the second is evidence.
	 */
	private boolean inValidatedPeak = false;
	/** Mass difference after subtracting the inferred isotope offset. */
	private double normalizedMassDiff = Double.NaN;
	private double peakCentroid = Double.NaN;
	private double distanceToPeakCentroid = Double.NaN;
	private double distanceInSigmas = Double.NaN;
	private double localDensity;
	private int peakIsotopeOffset;
	private double peakSigma = Double.NaN;
	private double peakAmplitude = Double.NaN;
	private double peakR2 = Double.NaN;

	// --- aggregate counts (computed once from raw data, never from scores) --
	private int samePeptideCount;
	private int sameProteinCount;
	private int samePeakCount;

	/**
	 * Features computed outside this class - MSBooster's predicted-spectrum and
	 * retention-time columns, and stage-1's derived ratios.
	 * <p>
	 * Held as a map rather than as fields so a feature set can grow without
	 * changing {@link OpenSearchFeatures#FEATURE_NAMES}, whose length and order
	 * stage 2's vector depends on. {@code OpenSearchFeatures.toInstance} falls
	 * back to this map for any requested name it does not compute itself.
	 */
	private java.util.Map<String, Double> externalFeatures;

	// --- validation output -------------------------------------------------
	private String engineSupport;
	private double classScore = Double.NaN;
	private double classLogit = Double.NaN;
	private double calibratedScore = Double.NaN;
	private double qValue = Double.NaN;

	/**
	 * Evidence that the delta mass is a real modification rather than a
	 * precursor mass error: how much better the spectrum is explained once
	 * fragment ions are shifted by the delta mass. This is the quantity the
	 * original implementation reconstructed by hand from raw spectra.
	 *
	 * @return NaN when MSFragger did not report both scores
	 */
	public double getDeltaMassScoreGain() {
		if (Double.isNaN(scoreWithDeltaMass) || Double.isNaN(scoreWithoutDeltaMass)) {
			return Double.NaN;
		}
		return scoreWithDeltaMass - scoreWithoutDeltaMass;
	}

	/**
	 * Localization confidence: how far the best site beats the runner-up.
	 *
	 * @return NaN when MSFragger did not report both scores
	 */
	public double getLocalizationDeltaScore() {
		if (Double.isNaN(scoreWithDeltaMass) || Double.isNaN(secondBestScoreWithDeltaMass)) {
			return Double.NaN;
		}
		return scoreWithDeltaMass - secondBestScoreWithDeltaMass;
	}

	/**
	 * Number of candidate sites MSFragger marked (lowercased) in
	 * {@link #getLocalizedSequence()} — higher means more ambiguous placement.
	 *
	 * @return -1 when no localized sequence was reported
	 */
	public int getLocalizationSiteCount() {
		if (localizedSequence == null || localizedSequence.isEmpty()) {
			return -1;
		}
		int n = 0;
		for (int i = 0; i < localizedSequence.length(); i++) {
			if (Character.isLowerCase(localizedSequence.charAt(i))) {
				n++;
			}
		}
		return n;
	}

	/**
	 * Cached, group-resolved value of {@link #getCanonicalModifiedSequence()}.
	 * Set exactly once, for every PSM sharing one {@link #getModificationFamilyKey()},
	 * by {@code OpenSearchPostAnalysisTsvExporter.resolveCanonicalModifiedSequences}
	 * - see that method's javadoc for why this must be resolved once per group
	 * rather than computed independently per PSM.
	 */
	private String resolvedModifiedSequence;

	/**
	 * The one canonical string that identifies "this peptide form" everywhere
	 * this PSM's evidence is joined against peptide-level data: the modified
	 * sequence carrying the delta mass, e.g. {@code PEPT[+14.0157]IDE}.
	 * <p>
	 * This is the single source of truth for that join key. Before this
	 * method existed, {@code OpenSearchMetaPsm} built its own {@code MetaPSM}
	 * full-sequence from {@link #getLocalizedSequence()} directly - MSFragger's
	 * lowercase-marked sequence, with no delta-mass suffix and no
	 * isotope-family collapsing - while {@code OpenSearchMetaPeptide} keyed
	 * itself on this exact canonical form. The two strings coincide for an
	 * unmodified PSM (both reduce to the bare peptide) but diverge for a real
	 * modification, so every PSM-to-peptide join keyed on {@code getModSeq()}
	 * - {@code pepPsmCountMap}, {@code pepExpectMap}, individual PEP,
	 * {@code DdaMagGenomeScoringService}'s per-protein intensity weighting -
	 * silently missed modified peptides. It went unnoticed because Phase 1 set
	 * every peptide's intensity to zero, so a missed intensity lookup and a
	 * successful one looked identical; wiring in real FlashLFQ quantification
	 * would have made the miss silently zero out intensity weighting for every
	 * modified peptide. Both sides now compute this one string.
	 * <p>
	 * <b>The string must never embed a per-spectrum localization guess as
	 * part of a peptide's identity.</b> An earlier version of this method used
	 * {@link #getLocalizedSequence()} - MSFragger's lowercase-marked sequence -
	 * as the string's base directly. Because that localization is a single
	 * spectrum's best guess at which residue carries the delta mass, and the
	 * winning site can flip between candidate residues when the score margin
	 * is thin, two spectra of the identical peptide backbone and identical
	 * delta mass could produce two different lowercase patterns - fragmenting
	 * one true peptide into multiple {@code validated_psms.tsv} rows and,
	 * consequently, multiple rows in FlashLFQ's {@code QuantifiedPeptides.tsv}.
	 * Localization ambiguity does not mean two chromatographically distinct
	 * species exist, so it must never change which rows a peptide's evidence
	 * is split across.
	 * <p>
	 * Returns the cached value {@code resolveCanonicalModifiedSequences} set
	 * for this PSM's group, which places the {@code [+mass]} bracket at the
	 * single representative (highest-scoring) PSM's localized site - decided
	 * once per group, so every member of the group reports the identical
	 * string regardless of its own individual localization guess. Falls back
	 * to {@link #getModificationFamilyKey()} - the same string with no
	 * position claim at all - when that resolution has not run; this fallback
	 * is always safe (unique, uppercase-only outside its brackets) even
	 * though it omits the representative's chosen position.
	 *
	 * @see uk.ac.quadram.metapilot.core.runner.OpenSearchPostAnalysisTsvExporter
	 */
	public String getCanonicalModifiedSequence() {
		return resolvedModifiedSequence != null ? resolvedModifiedSequence : getModificationFamilyKey();
	}

	/**
	 * The position-agnostic half of {@link #getCanonicalModifiedSequence()}:
	 * peptide backbone plus which delta-mass family this PSM belongs to, but
	 * never a specific localized residue. This is what grouping must use -
	 * see {@link #getCanonicalModifiedSequence()}'s javadoc for why keying on
	 * {@link #getLocalizedSequence()} instead silently over-fragmented
	 * peptides. Uppercase peptide letters only; any annotation lives inside
	 * the {@code []}/{@code [unexplained:...]} tag.
	 */
	public String getModificationFamilyKey() {
		if (isShiftLocalisedExplained()) {
			// Identification v2: the fragments placed the shift and an explanation fits there. The
			// family is the explanation's exact mass, so every PSM of this peptide and explanation
			// shares one row, whatever its observed mass or peak - and whichever search found it.
			return peptide + String.format(Locale.ROOT, "[%+.4f]", shiftLocalisation.explainedMass);
		}
		if (closedStage) {
			return getClosedModifiedSequence();
		}
		if (isShiftFoldedUnmodified()) {
			// An accepted-peak shift the fragments do not carry: the unmodified peptide (design doc, section 8).
			return peptide;
		}
		if (peakId < 0) {
			// No fitted hypothesis exists for this mass. Keep it explicitly
			// non-canonical until the reporting-tier decision is implemented; it
			// must not masquerade as an unmodified peptide.
			return peptide + String.format(Locale.ROOT, "[unexplained:%+.4f]", massDiff);
		}

		double canonicalMass = peakCentroid;
		if (Double.isNaN(canonicalMass)) {
			canonicalMass = getNormalizedMassDiff();
		}
		if (Math.abs(canonicalMass) <= 0.05) {
			// Isotope-triggered members of the zero family are the same unmodified
			// analyte. Do not preserve lowercase localization emitted for an
			// isotope-selection artefact.
			return peptide;
		}
		return peptide + String.format(Locale.ROOT, "[%+.4f]", canonicalMass);
	}

	public String getResolvedModifiedSequence() {
		return resolvedModifiedSequence;
	}

	/**
	 * Set once per {@link #getModificationFamilyKey()} group by
	 * {@code OpenSearchPostAnalysisTsvExporter.resolveCanonicalModifiedSequences},
	 * never per-PSM independently - that is what {@link #getCanonicalModifiedSequence()}
	 * now returns.
	 */
	public void setResolvedModifiedSequence(String resolvedModifiedSequence) {
		this.resolvedModifiedSequence = resolvedModifiedSequence;
	}

	/**
	 * Index (0-based) of the first residue {@link #getLocalizedSequence()}
	 * marks lowercase, or -1 when there is no localized sequence or no marked
	 * site. Used only to decide where the representative PSM of a group's
	 * {@code [+mass]} bracket goes - see
	 * {@code OpenSearchPostAnalysisTsvExporter#resolveCanonicalModifiedSequences}.
	 */
	public int getFirstLocalizedSiteIndex() {
		if (localizedSequence == null) {
			return -1;
		}
		for (int i = 0; i < localizedSequence.length(); i++) {
			if (Character.isLowerCase(localizedSequence.charAt(i))) {
				return i;
			}
		}
		return -1;
	}

	/**
	 * Cached, group-resolved value of {@link #getFlashLfqMonoisotopicMass()}.
	 * Set exactly once per {@link #getModificationFamilyKey()} group, by the
	 * same {@code OpenSearchPostAnalysisTsvExporter.resolveCanonicalModifiedSequences}
	 * pass that resolves {@link #getCanonicalModifiedSequence()} - it must be
	 * the same representative PSM and the same delta mass feeding both, or the
	 * numeric mass written to {@code validated_psms.tsv} would disagree with
	 * the mass the {@code [+mass]} bracket in "Full Sequence" implies.
	 * {@link Double#NaN} until resolved.
	 */
	private double resolvedMonoisotopicMass = Double.NaN;

	/**
	 * 0-based residue carrying the delta mass in {@link #getCanonicalModifiedSequence()},
	 * or -1 when the canonical form places it nowhere (unexplained, isotope-zero, or a
	 * representative without localization). Needed because the canonical string alone
	 * is ambiguous: {@code PEPTIDE[+x]} is both "unlocalized" (the family key) and
	 * "localized on the last residue". Set with the canonical sequence;
	 * {@link Integer#MIN_VALUE} until then.
	 */
	private int resolvedSiteIndex = Integer.MIN_VALUE;

	public int getResolvedSiteIndex() {
		return resolvedSiteIndex;
	}

	public void setResolvedSiteIndex(int resolvedSiteIndex) {
		this.resolvedSiteIndex = resolvedSiteIndex;
	}

	/**
	 * The mass to report to FlashLFQ as "Peptide Monoisotopic Mass" - the
	 * value it treats as this identification's theoretical precursor mass and
	 * searches the MS1 spectra for, at the stated retention time and charge.
	 * <p>
	 * {@link #getCalcNeutralPepMass()} is deliberately <b>not</b> this value.
	 * It is MSFragger's theoretical mass for the searched peptide string with
	 * no open delta applied - {@link #getMassDiff()} is defined as
	 * {@code precursor_neutral_mass - calc_neutral_pep_mass}, which is the
	 * entire mechanism of an open search: the engine deliberately does not
	 * fold the modification into the mass it calculates, because finding that
	 * residual is the point. Reporting the bare {@code calc_neutral_pep_mass}
	 * to FlashLFQ sends it looking for the precursor ion exactly
	 * {@code massDiff} Da away from where it actually sits in the MS1 data -
	 * for a real +79.9663 Da phosphorylation, over 10,000x FlashLFQ's own
	 * {@code PpmTolerance}. This is not a hypothetical: on the 12mix
	 * {@code combined_msfragger_sage} run, 299,591 of 896,297 validated PSM
	 * rows (33%) had {@code |massDiff| > 0.5 Da}, and every one of them would
	 * have been searched for at the wrong mass. This is a separate defect from
	 * the retention-time unit bug {@code retentionTimeMinutes} fixes - the two
	 * compound, but either alone accounts for FlashLFQ finding almost nothing.
	 * <p>
	 * Was, until this fix, copied unchanged from pFind's convention (see
	 * {@code MagOpenSearchHeadlessRunner.writeValidatedPsms}'s reference to
	 * {@code OpenPsmTsvWriter.addPFindPSM}), where it is correct: a closed
	 * search's calculated mass already includes whatever fixed/variable
	 * modifications the engine applied, so there is no undiscovered delta left
	 * to add. That precondition does not hold for an open search, where the
	 * whole delta mass is by definition excluded from
	 * {@code calc_neutral_pep_mass}.
	 * <p>
	 * Returns the cached value the resolution pass set for this PSM's group -
	 * {@link #getCalcNeutralPepMass()} plus that group's canonical delta mass,
	 * the same delta {@link #getCanonicalModifiedSequence()}'s bracket shows -
	 * so every PSM sharing one "Full Sequence" string also reports the
	 * identical mass, which FlashLFQ's own reader requires (it silently drops
	 * any later row whose mass does not exactly match the first one it saw for
	 * that Full Sequence). Falls back to this PSM's own calc_neutral_pep_mass
	 * plus its own delta when resolution has not run.
	 *
	 * @see uk.ac.quadram.metapilot.core.runner.OpenSearchPostAnalysisTsvExporter
	 */
	public double getFlashLfqMonoisotopicMass() {
		if (!Double.isNaN(resolvedMonoisotopicMass)) {
			return resolvedMonoisotopicMass;
		}
		if (isShiftLocalisedExplained()) {
			return calcNeutralPepMass + shiftLocalisation.explainedMass;
		}
		if (closedStage) {
			// A closed search's calculated mass already carries its variable modifications.
			return calcNeutralPepMass;
		}
		if (isShiftFoldedUnmodified()) {
			return calcNeutralPepMass;
		}
		if (peakId < 0) {
			return precursorNeutralMass;
		}
		double canonicalMass = peakCentroid;
		if (Double.isNaN(canonicalMass)) {
			canonicalMass = getNormalizedMassDiff();
		}
		if (Math.abs(canonicalMass) <= 0.05) {
			return calcNeutralPepMass;
		}
		return calcNeutralPepMass + canonicalMass;
	}

	public double getResolvedMonoisotopicMass() {
		return resolvedMonoisotopicMass;
	}

	/**
	 * Set once per {@link #getModificationFamilyKey()} group by
	 * {@code OpenSearchPostAnalysisTsvExporter.resolveCanonicalModifiedSequences},
	 * alongside {@link #setResolvedModifiedSequence(String)} and from the same
	 * representative PSM.
	 */
	public void setResolvedMonoisotopicMass(double resolvedMonoisotopicMass) {
		this.resolvedMonoisotopicMass = resolvedMonoisotopicMass;
	}

	public double getMassDiffPpm() {
		if (calcNeutralPepMass <= 0) {
			return Double.NaN;
		}
		return massDiff / calcNeutralPepMass * 1e6;
	}

	public double getNormalizedMassDiffPpm() {
		if (calcNeutralPepMass <= 0) {
			return Double.NaN;
		}
		return getNormalizedMassDiff() / calcNeutralPepMass * 1e6;
	}

	public int getPeptideLength() {
		return peptide == null ? 0 : peptide.length();
	}

	public double getMatchedIonFraction() {
		if (totNumIons <= 0) {
			return Double.NaN;
		}
		return (double) numMatchedIons / (double) totNumIons;
	}

	/** True when the fitter placed this PSM inside a delta-mass peak. */
	public boolean isInPeak() {
		return peakId >= 0;
	}

	/** True when that peak also passed the held-out binomial gate. */
	public boolean isInValidatedPeak() {
		return inValidatedPeak;
	}

	public void setInValidatedPeak(boolean inValidatedPeak) {
		this.inValidatedPeak = inValidatedPeak;
	}

	// --- accessors ---------------------------------------------------------

	public String getFileName() {
		return fileName;
	}

	public void setFileName(String fileName) {
		this.fileName = fileName;
	}

	public int getScan() {
		return scan;
	}

	public void setScan(int scan) {
		this.scan = scan;
	}

	public int getCharge() {
		return charge;
	}

	public void setCharge(int charge) {
		this.charge = charge;
	}

	public double getPrecursorNeutralMass() {
		return precursorNeutralMass;
	}

	public void setPrecursorNeutralMass(double precursorNeutralMass) {
		this.precursorNeutralMass = precursorNeutralMass;
	}

	/** Seconds, for both engines (SageOpenSearchReader converts Sage's minutes). */
	public double getRetentionTime() {
		return retentionTime;
	}

	public void setRetentionTime(double retentionTime) {
		this.retentionTime = retentionTime;
	}

	public int getHitRank() {
		return hitRank;
	}

	public void setHitRank(int hitRank) {
		this.hitRank = hitRank;
	}

	public String getPeptide() {
		return peptide;
	}

	public void setPeptide(String peptide) {
		this.peptide = peptide;
	}

	public String getPrevAA() {
		return prevAA;
	}

	public void setPrevAA(String prevAA) {
		this.prevAA = prevAA;
	}

	public String getNextAA() {
		return nextAA;
	}

	public void setNextAA(String nextAA) {
		this.nextAA = nextAA;
	}

	public double getCalcNeutralPepMass() {
		return calcNeutralPepMass;
	}

	public void setCalcNeutralPepMass(double calcNeutralPepMass) {
		this.calcNeutralPepMass = calcNeutralPepMass;
	}

	public double getMassDiff() {
		return massDiff;
	}

	public void setMassDiff(double massDiff) {
		this.massDiff = massDiff;
	}

	public int getNumTolTerm() {
		return numTolTerm;
	}

	public void setNumTolTerm(int numTolTerm) {
		this.numTolTerm = numTolTerm;
	}

	public int getNumMissedCleavages() {
		return numMissedCleavages;
	}

	public void setNumMissedCleavages(int numMissedCleavages) {
		this.numMissedCleavages = numMissedCleavages;
	}

	public int getNumMatchedIons() {
		return numMatchedIons;
	}

	public void setNumMatchedIons(int numMatchedIons) {
		this.numMatchedIons = numMatchedIons;
	}

	public int getTotNumIons() {
		return totNumIons;
	}

	public void setTotNumIons(int totNumIons) {
		this.totNumIons = totNumIons;
	}

	public String getProtein() {
		return protein;
	}

	public void setProtein(String protein) {
		this.protein = protein;
	}

	/** Delegate plus alternative proteins; never null, may be empty. */
	public String[] getAllProteins() {
		return allProteins;
	}

	public void setAllProteins(String[] allProteins) {
		this.allProteins = allProteins == null ? new String[0] : allProteins;
	}

	public boolean isTarget() {
		return target;
	}

	public void setTarget(boolean target) {
		this.target = target;
	}

	public double getHyperscore() {
		return hyperscore;
	}

	public void setHyperscore(double hyperscore) {
		this.hyperscore = hyperscore;
	}

	public double getNextscore() {
		return nextscore;
	}

	public void setNextscore(double nextscore) {
		this.nextscore = nextscore;
	}

	public double getExpect() {
		return expect;
	}

	public void setExpect(double expect) {
		this.expect = expect;
	}

	public String getLocalizedSequence() {
		return localizedSequence;
	}

	/**
	 * Why this PSM's delta mass was moved, or null when it was left as searched.
	 * <p>
	 * {@code PRECURSOR_ARTIFACT} - {@link DeltaMassEvidenceGate} found no
	 * shifted-fragment support for the peak this PSM sat in and folded it onto
	 * the zero family; the peptide is unmodified and the precursor mass was
	 * wrong. {@code SUBSTITUTION:<from><pos><to>} -
	 * {@link SubstitutionReassigner} found the delta is explained by a single
	 * amino-acid change at a position the localization agrees with, so the
	 * sequence rather than the mass is what is wrong.
	 * <p>
	 * {@link #getMassDiff()} is never rewritten by either: it stays the record of
	 * what the instrument reported, and only {@link #getNormalizedMassDiff()}
	 * moves.
	 */
	public String getDeltaMassReassignment() {
		return deltaMassReassignment;
	}

	public void setDeltaMassReassignment(String deltaMassReassignment) {
		this.deltaMassReassignment = deltaMassReassignment;
	}

	/**
	 * {@link SubstitutionReassigner.Verdict} as a name, or null when this PSM was
	 * not on a substitution-shaped peak.
	 */
	public String getSubstitutionVerdict() {
		return substitutionVerdict;
	}

	public void setSubstitutionVerdict(String substitutionVerdict) {
		this.substitutionVerdict = substitutionVerdict;
	}

	/**
	 * The sequence this spectrum probably came from, if a single-residue change
	 * explains its delta mass. A proposal, not a decision: confirming it needs the
	 * spectrum scored against it, which only the search engine can do.
	 */
	public String getSubstitutionCandidate() {
		return substitutionCandidate;
	}

	public void setSubstitutionCandidate(String substitutionCandidate) {
		this.substitutionCandidate = substitutionCandidate;
	}

	public void setLocalizedSequence(String localizedSequence) {
		this.localizedSequence = localizedSequence;
	}

	public double getScoreWithoutDeltaMass() {
		return scoreWithoutDeltaMass;
	}

	public void setScoreWithoutDeltaMass(double scoreWithoutDeltaMass) {
		this.scoreWithoutDeltaMass = scoreWithoutDeltaMass;
	}

	public double getScoreWithDeltaMass() {
		return scoreWithDeltaMass;
	}

	public void setScoreWithDeltaMass(double scoreWithDeltaMass) {
		this.scoreWithDeltaMass = scoreWithDeltaMass;
	}

	public double getSecondBestScoreWithDeltaMass() {
		return secondBestScoreWithDeltaMass;
	}

	public void setSecondBestScoreWithDeltaMass(double secondBestScoreWithDeltaMass) {
		this.secondBestScoreWithDeltaMass = secondBestScoreWithDeltaMass;
	}

	public double getIonsWithoutDeltaMass() {
		return ionsWithoutDeltaMass;
	}

	public void setIonsWithoutDeltaMass(double ionsWithoutDeltaMass) {
		this.ionsWithoutDeltaMass = ionsWithoutDeltaMass;
	}

	public double getIonsWithDeltaMass() {
		return ionsWithDeltaMass;
	}

	public void setIonsWithDeltaMass(double ionsWithDeltaMass) {
		this.ionsWithDeltaMass = ionsWithDeltaMass;
	}

	public Map<String, String> getSearchScores() {
		return searchScores;
	}

	public void setSearchScores(Map<String, String> searchScores) {
		this.searchScores = searchScores == null ? Collections.emptyMap() : new LinkedHashMap<>(searchScores);
	}

	public Map<String, String> getPtmResultAttributes() {
		return ptmResultAttributes;
	}

	public void setPtmResultAttributes(Map<String, String> ptmResultAttributes) {
		this.ptmResultAttributes = ptmResultAttributes == null ? Collections.emptyMap()
				: new LinkedHashMap<>(ptmResultAttributes);
	}

	public int getPeakId() {
		return peakId;
	}

	public void setPeakId(int peakId) {
		this.peakId = peakId;
	}

	public int getObservedPeakId() {
		return observedPeakId;
	}

	public void setObservedPeakId(int observedPeakId) {
		this.observedPeakId = observedPeakId;
	}

	public double getNormalizedMassDiff() {
		return Double.isNaN(normalizedMassDiff) ? massDiff : normalizedMassDiff;
	}

	public void setNormalizedMassDiff(double normalizedMassDiff) {
		this.normalizedMassDiff = normalizedMassDiff;
	}

	public double getPeakCentroid() {
		return peakCentroid;
	}

	public void setPeakCentroid(double peakCentroid) {
		this.peakCentroid = peakCentroid;
	}

	public double getDistanceToPeakCentroid() {
		return distanceToPeakCentroid;
	}

	public void setDistanceToPeakCentroid(double distanceToPeakCentroid) {
		this.distanceToPeakCentroid = distanceToPeakCentroid;
	}

	public double getDistanceInSigmas() {
		return distanceInSigmas;
	}

	public void setDistanceInSigmas(double distanceInSigmas) {
		this.distanceInSigmas = distanceInSigmas;
	}

	public double getLocalDensity() {
		return localDensity;
	}

	public void setLocalDensity(double localDensity) {
		this.localDensity = localDensity;
	}

	public int getPeakIsotopeOffset() {
		return peakIsotopeOffset;
	}

	public void setPeakIsotopeOffset(int peakIsotopeOffset) {
		this.peakIsotopeOffset = peakIsotopeOffset;
	}

	public double getPeakSigma() {
		return peakSigma;
	}

	public void setPeakSigma(double peakSigma) {
		this.peakSigma = peakSigma;
	}

	public double getPeakAmplitude() {
		return peakAmplitude;
	}

	public void setPeakAmplitude(double peakAmplitude) {
		this.peakAmplitude = peakAmplitude;
	}

	public double getPeakR2() {
		return peakR2;
	}

	public void setPeakR2(double peakR2) {
		this.peakR2 = peakR2;
	}

	public int getSamePeptideCount() {
		return samePeptideCount;
	}

	public void setSamePeptideCount(int samePeptideCount) {
		this.samePeptideCount = samePeptideCount;
	}

	public int getSameProteinCount() {
		return sameProteinCount;
	}

	public void setSameProteinCount(int sameProteinCount) {
		this.sameProteinCount = sameProteinCount;
	}

	/** The named external feature, or null when this PSM has no value for it. */
	public Double getExternalFeature(String name) {
		return externalFeatures == null || name == null ? null : externalFeatures.get(name);
	}

	public java.util.Map<String, Double> getExternalFeatures() {
		return externalFeatures == null ? java.util.Collections.emptyMap() : externalFeatures;
	}

	public void setExternalFeatures(java.util.Map<String, Double> externalFeatures) {
		this.externalFeatures = externalFeatures;
	}

	public int getSamePeakCount() {
		return samePeakCount;
	}

	public void setSamePeakCount(int samePeakCount) {
		this.samePeakCount = samePeakCount;
	}

	public double getClassScore() {
		return classScore;
	}

	/**
	 * Which engines proposed this candidate: MSFRAGGER, SAGE, or BOTH.
	 * <p>
	 * Distinct from {@link #getEngine()}, which names the engine that produced
	 * this particular row. When engine rows are merged into one combined
	 * candidate only one row is carried downstream, so {@code engine} alone would
	 * report every jointly-proposed identification as coming from whichever row
	 * happened to be primary - making "what did the second engine contribute"
	 * unanswerable from the output, which is the question the combined path
	 * exists to let you ask.
	 */
	public String getEngineSupport() {
		return engineSupport == null ? String.valueOf(engine) : engineSupport;
	}

	public void setEngineSupport(String engineSupport) {
		this.engineSupport = engineSupport;
	}

	public void setClassScore(double classScore) {
		this.classScore = classScore;
	}

	public double getClassLogit() {
		return classLogit;
	}

	public void setClassLogit(double classLogit) {
		this.classLogit = classLogit;
	}

	/**
	 * The engine-comparable score {@link PerEngineClassifier} assigns, NaN on the
	 * combined-candidate path. It saturates at the cap for most accepted PSMs, so
	 * it serves only to rank PSMs from different engines against each other.
	 */
	public double getCalibratedScore() {
		return calibratedScore;
	}

	public void setCalibratedScore(double calibratedScore) {
		this.calibratedScore = calibratedScore;
	}

	/**
	 * The score for ranking PSMs that may come from different engines: the
	 * calibrated score when set, otherwise the class logit, then the class score.
	 */
	public double getComparableScore() {
		if (!Double.isNaN(calibratedScore)) {
			return calibratedScore;
		}
		return Double.isNaN(classLogit) ? classScore : classLogit;
	}

	/**
	 * Joint closed+open identification only: the FDR group this PSM was accepted or walked in
	 * ({@code JointCandidate.Group} name), null on every other path. Written to
	 * {@code validated_psms.tsv} as {@code id_group}.
	 */
	private String identificationGroup;

	public String getIdentificationGroup() {
		return identificationGroup;
	}

	public void setIdentificationGroup(String identificationGroup) {
		this.identificationGroup = identificationGroup;
	}

	/**
	 * Identification v2 only: the {@link ShiftLocaliser} result of an unexplained-shift spectrum
	 * winner, null otherwise. Feeds the pangenome variant check's localised candidates.
	 */
	private ShiftLocaliser.Result shiftLocalisation;

	public ShiftLocaliser.Result getShiftLocalisation() {
		return shiftLocalisation;
	}

	public void setShiftLocalisation(ShiftLocaliser.Result shiftLocalisation) {
		this.shiftLocalisation = shiftLocalisation;
	}

	/**
	 * Whether this PSM's shift is explained - localised by the fragments and explained there, or an
	 * accepted peak's unplaced shift given its UNIMOD modification ({@code ShiftLocaliser.State#ASSIGNED}):
	 * its exported sequence then carries the explanation's mass on that residue, and FlashLFQ is given
	 * the peptide's mass plus it. Open PSMs, and closed PSMs whose isotope offset was
	 * an established shift ({@code JointIdentification.localiseWinners}); a closed PSM carrying variable
	 * modifications is never localised.
	 */
	public boolean isShiftLocalisedExplained() {
		return shiftLocalisation != null && shiftLocalisation.isExplained()
				&& !Double.isNaN(shiftLocalisation.explainedMass)
				&& (!closedStage || getClosedModificationCount() == 0);
	}

	/**
	 * An accepted-shift PSM (on an accepted peak) whose shift the fragments do not carry: exported as the
	 * unmodified peptide, the observed shift kept in its mass-difference columns.
	 */
	public boolean isShiftFoldedUnmodified() {
		return !closedStage && shiftLocalisation != null && shiftLocalisation.state == ShiftLocaliser.State.UNSUPPORTED
				&& "OPEN_ACCEPTED_SHIFT".equals(identificationGroup);
	}

	/**
	 * The site of an exported shift: the fragment-supported one (explained, or not explained there), or
	 * the UNIMOD site an unplaced accepted shift was assigned; -1 otherwise.
	 */
	public int getLocalisedShiftSite() {
		if (shiftLocalisation == null || shiftLocalisation.site < 0
				|| closedStage && getClosedModificationCount() != 0) {
			return -1;
		}
		return shiftLocalisation.state == ShiftLocaliser.State.LOCALISED_EXPLAINED
				|| shiftLocalisation.state == ShiftLocaliser.State.LOCALISED_UNEXPLAINED
				|| shiftLocalisation.state == ShiftLocaliser.State.ASSIGNED ? shiftLocalisation.site : -1;
	}

	/**
	 * The explanation of the representative PSM that placed this PSM's group bracket, set with
	 * {@link #getCanonicalModifiedSequence()}; empty for groups that are not localised-explained.
	 */
	private String resolvedShiftExplanation = "";

	public String getResolvedShiftExplanation() {
		return resolvedShiftExplanation;
	}

	public void setResolvedShiftExplanation(String resolvedShiftExplanation) {
		this.resolvedShiftExplanation = resolvedShiftExplanation == null ? "" : resolvedShiftExplanation;
	}

	/**
	 * Beside {@link #getResolvedShiftExplanation()}, for the peptide table: how the representative PSM's
	 * site is known ({@code localised_explained}: placed by the fragments; {@code assigned}: the UNIMOD
	 * site of the peak's mass) and what kind of explanation it is ({@code modification}, ...). Empty
	 * without an explanation.
	 */
	private String resolvedShiftState = "";
	private String resolvedShiftClass = "";

	public String getResolvedShiftState() {
		return resolvedShiftState;
	}

	public String getResolvedShiftClass() {
		return resolvedShiftClass;
	}

	public void setResolvedShiftStateAndClass(String state, String shiftClass) {
		this.resolvedShiftState = state == null ? "" : state;
		this.resolvedShiftClass = shiftClass == null ? "" : shiftClass;
	}

	public double getQValue() {
		return qValue;
	}

	public void setQValue(double qValue) {
		this.qValue = qValue;
	}

	@Override
	public String toString() {
		return fileName + ":" + scan + " rank=" + hitRank + " " + peptide + " dm=" + massDiff
				+ (target ? " target" : " decoy");
	}
}
