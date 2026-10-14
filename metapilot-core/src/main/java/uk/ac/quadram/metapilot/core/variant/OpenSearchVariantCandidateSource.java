package uk.ac.quadram.metapilot.core.variant;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.IdentityHashMap;
import java.util.List;
import java.util.Locale;
import java.util.Set;
import java.util.TreeSet;

import uk.ac.quadram.metapilot.core.opensearch.IsobaricModificationTable;
import uk.ac.quadram.metapilot.core.opensearch.OpenSearchPsm;
import uk.ac.quadram.metapilot.core.opensearch.SubstitutionReassigner;

/**
 * Front end for the MetaPilot open search: substitution candidates from what
 * {@link SubstitutionReassigner} already proposed.
 * <p>
 * Two sources, kept apart by {@code peak_substitution_shaped}:
 * <ul>
 * <li>{@code true} - proposals on peaks the reassigner judged
 * substitution-shaped; these PSMs carry its verdict.</li>
 * <li>{@code false} - proposals it computed on peaks whose centroid matches a
 * substitution mass but which failed its shape or in-database test
 * ({@link SubstitutionReassigner.Result#rejectedPeakProposals}). Nothing in the
 * workflow acts on them; they are here so the modification-like peaks
 * (+0.9841, +40.0313, +57.0214, ...) can be compared with the
 * substitution-shaped ones on the same footing.</li>
 * </ul>
 * Only proposals on validated target PSMs become candidates: the reassigner
 * runs on every candidate before FDR, and the check is about identifications
 * the workflow reports.
 * <p>
 * Masses are the reassigner's own ({@code SubstitutionMassTable}, cysteine
 * carbamidomethylated as the search fixed it); MS2 decidability uses the
 * shared {@link IsobaricModificationTable} with
 * {@code substitution_mass_tolerance_da}, the tolerance that already feeds
 * {@code sequence_variants.tsv}'s isobaric column.
 */
public final class OpenSearchVariantCandidateSource {

	public static final class Result {
		public final List<VariantCandidate> candidates = new ArrayList<>();
		/** Peptides of every validated target PSM; the matched-control pool. */
		public final Set<String> basePeptides = new TreeSet<>();
		public int shapedCandidates;
		public int rejectedPeakCandidates;

		public String describe() {
			return String.format(Locale.ROOT, "%,d candidates (%,d on substitution-shaped peaks, %,d on rejected "
					+ "peaks) from %,d validated target peptides", candidates.size(), shapedCandidates,
					rejectedPeakCandidates, basePeptides.size());
		}
	}

	private OpenSearchVariantCandidateSource() {
	}

	/**
	 * @param substitutions the reassigner's result for this validation; null
	 *                      gives no candidates
	 * @param validated     the validated PSMs, compared by identity
	 */
	public static Result collect(SubstitutionReassigner.Result substitutions, Collection<OpenSearchPsm> validated,
			double undecidableToleranceDa) {
		Result result = new Result();
		Set<OpenSearchPsm> accepted = Collections.newSetFromMap(new IdentityHashMap<>());
		for (OpenSearchPsm psm : validated) {
			if (psm.isTarget() && psm.getPeptide() != null && !psm.getPeptide().isEmpty()) {
				accepted.add(psm);
				result.basePeptides.add(psm.getPeptide());
			}
		}
		if (substitutions == null) {
			return result;
		}
		for (SubstitutionReassigner.Proposal p : substitutions.proposals) {
			if (accepted.contains(p.psm)) {
				result.candidates.add(candidate(p, true, undecidableToleranceDa));
				result.shapedCandidates++;
			}
		}
		for (SubstitutionReassigner.Proposal p : substitutions.rejectedPeakProposals) {
			if (accepted.contains(p.psm)) {
				result.candidates.add(candidate(p, false, undecidableToleranceDa));
				result.rejectedPeakCandidates++;
			}
		}
		return result;
	}

	static VariantCandidate candidate(SubstitutionReassigner.Proposal p, boolean peakShaped,
			double undecidableToleranceDa) {
		OpenSearchPsm psm = p.psm;
		VariantCandidate c = new VariantCandidate();
		c.engine = VariantCandidate.Engine.OPEN_SEARCH;
		c.file = psm.getFileName();
		c.scan = psm.getScan();
		c.psmId = psm.getHitRank();
		c.charge = psm.getCharge();
		c.qValue = psm.getQValue();
		if (c.charge > 0 && psm.getPrecursorNeutralMass() > 0) {
			c.precursorMz = psm.getPrecursorNeutralMass() / c.charge + PROTON;
		}
		for (String protein : psm.getAllProteins()) {
			if (protein != null && !protein.isEmpty()) {
				c.assignedProteins.add(protein);
			}
		}
		if (c.assignedProteins.isEmpty() && psm.getProtein() != null && !psm.getProtein().isEmpty()) {
			c.assignedProteins.add(psm.getProtein());
		}
		c.basePeptide = psm.getPeptide();
		c.position = p.position + 1;
		c.fromAa = p.substitution.from;
		c.toAa = p.substitution.to;
		c.variantPeptide = p.candidate;
		c.deltaMass = p.substitution.massDelta;
		c.peakId = psm.getPeakId();
		c.peakSubstitutionShaped = peakShaped;
		c.isobaricModifications.addAll(IsobaricModificationTable.sameResidue(c.fromAa, c.deltaMass, c.position == 1,
				undecidableToleranceDa));
		c.ms2Decidable = c.isobaricModifications.isEmpty();
		c.frontEndEvidence = peakShaped ? String.valueOf(psm.getSubstitutionVerdict()) : "peak_not_substitution_shaped";
		c.localizedSequence = psm.getLocalizedSequence() == null ? "" : psm.getLocalizedSequence();
		c.reassignerInDb = p.inDatabase;
		return c;
	}

	static final double PROTON = 1.007276466812;
}
