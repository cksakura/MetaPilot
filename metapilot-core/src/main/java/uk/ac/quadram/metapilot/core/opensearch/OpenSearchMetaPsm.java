package uk.ac.quadram.metapilot.core.opensearch;

import java.util.ArrayList;
import java.util.Collections;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

import uk.ac.quadram.metapilot.task.io.pep.MetaPSM;
import uk.ac.quadram.metapilot.task.io.qpx.EngineModificationParser;
import uk.ac.quadram.metapilot.task.io.qpx.QpxModification;
import uk.ac.quadram.metapilot.task.io.qpx.QpxScore;

/**
 * Adapts an {@link OpenSearchPsm} to {@link MetaPSM}, the shape the shared MAG
 * services consume.
 * <p>
 * This is the seam between the engine-specific segment of this workflow and the
 * generic MAG skeleton: once PSMs are in this form, genome selection
 * ({@code MagGenomeSelectionService}), genome scoring, protein inference and
 * the downstream taxonomy/function/report steps all work unchanged.
 * <p>
 * {@code score} carries the classifier logit where one has been assigned,
 * falling back to the classifier probability, then to hyperscore for stage-1
 * PSMs, which reach genome selection before any classifier has run.
 * <p>
 * The logit is preferred over the probability deliberately: {@code classScore}
 * is a sigmoid-compressed probability that saturates at/near 1.0 for the
 * large majority of confidently-accepted targets (measured: ~80% of targets
 * in a real run have {@code classScore >= 0.9999}), which destroys the
 * dynamic range that {@code DdaMagGenomeScoringService}'s top-k bootstrap
 * needs to separate real genome signal from decoy noise - and, worse, the
 * decoys that survive PSM-level FDR filtering are by construction the ones
 * whose probability sits right next to the accepted targets', so decoy and
 * target scores collapse together exactly where discrimination matters most.
 * {@code classLogit} (pre-sigmoid) retains that separation. The per-engine
 * calibrated score ({@link OpenSearchPsm#getCalibratedScore()}) must not be used
 * here: it is capped at 12 and most accepted PSMs reach the cap, which on the
 * 2026-10-08 IPX0003851001 MSFragger-only run left 2 of 32 genomes at
 * q <= 0.01 (31 with the logit; all 32 once the validated decoys are the
 * genome-scoring background). {@code pep} carries the q-value,
 * matching {@link MetaPSM#getqValue()}'s contract.
 */
public class OpenSearchMetaPsm extends MetaPSM {

	private final OpenSearchPsm source;

	/**
	 * Fixed modifications of the search, keyed the way
	 * {@code SageOpenSearchRunner.Settings.staticMods} keys them (a residue letter,
	 * {@code ^} peptide N-term, {@code $} peptide C-term). Neither engine's peptide
	 * string carries them, but calc_neutral_pep_mass includes them, so the QPX
	 * peptidoform must state them. Empty when unknown.
	 */
	private final Map<String, Double> staticMods;

	public OpenSearchMetaPsm(OpenSearchPsm source) {
		this(source, Collections.emptyMap());
	}

	public OpenSearchMetaPsm(OpenSearchPsm source, Map<String, Double> staticMods) {
		// The full/mod sequence is the join key every peptide-level lookup in
		// DdaMagPostAnalysisService and DdaMagGenomeScoringService uses
		// (pepPsmCountMap, pepExpectMap, modSeqToPeptide...), so it must be
		// built exactly the way OpenSearchMetaPeptide builds its own modSeq -
		// see OpenSearchPsm.getCanonicalModifiedSequence() for why this used
		// to be getLocalizedSequence() and why that silently broke the join
		// for every modified peptide.
		super(source.getFileName(), source.getScan(), source.getPeptide(), source.getCanonicalModifiedSequence(),
				source.getCharge(), proteinsOf(source), scoreOf(source),
				Double.isNaN(source.getQValue()) ? 1.0 : source.getQValue());
		this.source = source;
		this.staticMods = staticMods == null ? Collections.emptyMap() : staticMods;
	}

	private static double scoreOf(OpenSearchPsm psm) {
		// Prefer the pre-sigmoid logit, which keeps the dynamic range the genome-scoring
		// bootstrap needs; classScore alone saturates near 1.0 for most targets
		// and for the decoys that survive FDR filtering, collapsing exactly the
		// separation a top-k/bootstrap comparison depends on.
		if (!Double.isNaN(psm.getClassLogit())) {
			return psm.getClassLogit();
		}
		if (!Double.isNaN(psm.getClassScore())) {
			return psm.getClassScore();
		}
		// Stage-1 PSMs reach genome selection before the classifier exists.
		return Double.isNaN(psm.getHyperscore()) ? 0.0 : psm.getHyperscore();
	}

	/**
	 * The full protein list, matching what {@code FragpipePSMTsvReader} builds
	 * from "Protein" + "Mapped Proteins". Genome selection counts genomes per
	 * peptide, so passing the delegate alone would credit a shared marker peptide
	 * to a single genome and change which genomes survive the knee cut.
	 */
	private static String[] proteinsOf(OpenSearchPsm psm) {
		String[] all = psm.getAllProteins();
		if (all != null && all.length > 0) {
			return all;
		}
		String protein = psm.getProtein();
		if (protein == null || protein.isBlank()) {
			return new String[0];
		}
		return new String[] { protein.trim() };
	}

	public OpenSearchPsm getSource() {
		return source;
	}

	/** True for decoy hits; genome selection filters these out itself. */
	public boolean isDecoy() {
		return !source.isTarget();
	}

	@Override
	public double getExpMass() {
		return source.getPrecursorNeutralMass();
	}

	@Override
	public double getRt() {
		return source.getRetentionTime();
	}

	@Override
	public double getCalMass() {
		return source.getCalcNeutralPepMass();
	}

	@Override
	public double getMassDiff() {
		return source.getNormalizedMassDiff();
	}

	@Override
	public int getMiss() {
		return source.getNumMissedCleavages();
	}

	// ------------------------------------------------------------------
	// Normalised accessors (see MetaPSM)
	// ------------------------------------------------------------------

	@Override
	public String getSearchEngine() {
		return source.getEngine() == OpenSearchPsm.SearchEngine.SAGE ? "Sage" : "MSFragger";
	}

	@Override
	public double getRtSeconds() {
		return source.getRetentionTime();
	}

	@Override
	public double getObservedNeutralMass() {
		return source.getPrecursorNeutralMass();
	}

	/**
	 * calc_neutral_pep_mass plus the group's canonical delta - the mass of exactly the
	 * peptidoform {@link #getProForma()} writes (and what FlashLFQ is given). For an
	 * unexplained delta that is the observed mass, the delta being written unlocalized.
	 */
	@Override
	public double getTheoreticalNeutralMass() {
		return source.getFlashLfqMonoisotopicMass();
	}

	@Override
	public double getPsmQValue() {
		return source.getQValue();
	}

	/**
	 * Fixed modifications on their residues (mapped to UNIMOD by mass), plus the open
	 * delta mass as the canonical form states it: on the representative's residue when
	 * localized, otherwise unlocalized. Open deltas stay mass-only on purpose - naming
	 * them is the modification-verdict step's job, not a mass lookup's.
	 */
	@Override
	public List<QpxModification> getQpxModifications() {
		String seq = getBaseSequence();
		if (seq == null) {
			return null;
		}
		int n = seq.length();
		List<QpxModification> mods = new ArrayList<>();
		Double nTerm = staticMods.get("^");
		if (nTerm != null && nTerm != 0.0) {
			mods.add(EngineModificationParser.massMod(nTerm, seq, 0));
		}
		for (int i = 0; i < n; i++) {
			Double m = staticMods.get(String.valueOf(seq.charAt(i)));
			if (m != null && m != 0.0) {
				mods.add(EngineModificationParser.massMod(m, seq, i + 1));
			}
		}
		Double cTerm = staticMods.get("$");
		if (cTerm != null && cTerm != 0.0) {
			mods.add(EngineModificationParser.massMod(cTerm, seq, n + 1));
		}

		if (source.isClosedStage()) {
			double[] variable = source.getClosedModifications();
			for (int i = 0; variable != null && i < variable.length; i++) {
				if (variable[i] != 0.0) {
					mods.add(EngineModificationParser.massMod(variable[i], seq, i));
				}
			}
			return mods;
		}
		if (source.getPeakId() < 0) {
			double delta = source.getMassDiff();
			if (!Double.isNaN(delta)) {
				mods.add(QpxModification.unlocalized(null, delta, seq, null));
			}
			return mods;
		}
		double delta = getTheoreticalNeutralMass() - source.getCalcNeutralPepMass();
		if (Double.isNaN(delta) || Math.abs(delta) <= 0.05) {
			return mods;
		}
		int site = source.getResolvedSiteIndex();
		if (site >= 0 && site < n) {
			mods.add(QpxModification.at(null, delta, seq, site + 1));
		} else {
			mods.add(QpxModification.unlocalized(null, delta, seq, null));
		}
		return mods;
	}

	@Override
	public List<QpxScore> getAdditionalScores() {
		List<QpxScore> scores = new ArrayList<>(6);
		String prefix = source.getEngine() == OpenSearchPsm.SearchEngine.SAGE ? "sage" : "msfragger";
		QpxScore.addIfFinite(scores, prefix + "_hyperscore", source.getHyperscore(), Boolean.TRUE);
		QpxScore.addIfFinite(scores, prefix + "_expect", source.getExpect(), Boolean.FALSE);
		QpxScore.addIfFinite(scores, "metapilot_class_score", source.getClassScore(), Boolean.TRUE);
		QpxScore.addIfFinite(scores, "metapilot_class_logit", source.getClassLogit(), Boolean.TRUE);
		QpxScore.addIfFinite(scores, "q_value", source.getQValue(), Boolean.FALSE);
		QpxScore.addIfFinite(scores, "delta_mass", source.getMassDiff(), null);
		return scores;
	}

	@Override
	public Map<String, String> getCvParams() {
		Map<String, String> cv = new LinkedHashMap<>();
		cv.put("metapilot:delta_mass_peak_id", String.valueOf(source.getPeakId()));
		if (source.isClosedStage()) {
			cv.put("metapilot:identification_stage", "closed");
		} else if (source.getPeakId() < 0) {
			cv.put("metapilot:delta_mass_unexplained", "true");
		}
		putIfSet(cv, "metapilot:engine_support", source.getEngineSupport());
		putIfSet(cv, "metapilot:substitution_verdict", source.getSubstitutionVerdict());
		putIfSet(cv, "metapilot:substitution_candidate", source.getSubstitutionCandidate());
		return cv;
	}

	private static void putIfSet(Map<String, String> cv, String key, String value) {
		if (value != null && !value.isBlank()) {
			cv.put(key, value);
		}
	}
}
