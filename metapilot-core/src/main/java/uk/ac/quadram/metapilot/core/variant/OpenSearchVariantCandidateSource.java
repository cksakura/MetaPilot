package uk.ac.quadram.metapilot.core.variant;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.IdentityHashMap;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.Set;
import java.util.TreeSet;

import uk.ac.quadram.metapilot.core.opensearch.IsobaricModificationTable;
import uk.ac.quadram.metapilot.core.opensearch.OpenSearchPsm;
import uk.ac.quadram.metapilot.core.opensearch.ShiftLocaliser;
import uk.ac.quadram.metapilot.core.opensearch.SubstitutionMassTable;
import uk.ac.quadram.metapilot.core.opensearch.SubstitutionReassigner;

/**
 * Front end for the MetaPilot open search: substitution candidates from what
 * {@link SubstitutionReassigner} already proposed, and from shift localisation.
 * <p>
 * Three sources, kept apart by {@code peak_substitution_shaped}:
 * <ul>
 * <li>{@code true} - proposals on peaks the reassigner judged
 * substitution-shaped; these PSMs carry its verdict.</li>
 * <li>{@code false} - proposals it computed on peaks whose centroid matches a
 * substitution mass but which failed its shape or in-database test
 * ({@link SubstitutionReassigner.Result#rejectedPeakProposals}). Nothing in the
 * workflow acts on them; they are here so the modification-like peaks
 * (+0.9841, +40.0313, +57.0214, ...) can be compared with the
 * substitution-shaped ones on the same footing.</li>
 * <li>NA - identification v2: {@code OPEN_LOCALISED} and {@code OPEN_SUBSTITUTION_SITE} PSMs whose
 * {@link ShiftLocaliser} explanation at the localised site includes a
 * substitution, on a peak or not. Each such substitution is a candidate at
 * that site, with the site probability and fragment support;
 * {@code front_end_call} is {@value #LOCALISED_CALL}. Where a peak proposal
 * names the same PSM, site and residue, the two are one candidate carrying
 * both.</li>
 * </ul>
 * Only validated target PSMs give candidates: the reassigner runs on every
 * candidate before FDR, and the check is about identifications the workflow
 * reports.
 * <p>
 * Masses are the reassigner's own ({@code SubstitutionMassTable}, cysteine
 * carbamidomethylated as the search fixed it); MS2 decidability uses the
 * shared {@link IsobaricModificationTable} with
 * {@code substitution_mass_tolerance_da}, the tolerance that already feeds
 * {@code sequence_variants.tsv}'s isobaric column.
 */
public final class OpenSearchVariantCandidateSource {

	/** {@code front_end_call} of a candidate from shift localisation alone. */
	public static final String LOCALISED_CALL = "localised";
	private static final String SUBSTITUTION_LABEL = "subst ";

	public static final class Result {
		public final List<VariantCandidate> candidates = new ArrayList<>();
		/** Peptides of every validated target PSM; the matched-control pool. */
		public final Set<String> basePeptides = new TreeSet<>();
		public int shapedCandidates;
		public int rejectedPeakCandidates;
		/** Candidates from shift localisation alone; peak proposals it confirmed are not counted again. */
		public int localisedCandidates;

		public String describe() {
			return String.format(Locale.ROOT, "%,d candidates (%,d on substitution-shaped peaks, %,d on rejected "
					+ "peaks, %,d from shift localisation alone) from %,d validated target peptides",
					candidates.size(), shapedCandidates, rejectedPeakCandidates, localisedCandidates,
					basePeptides.size());
		}
	}

	private OpenSearchVariantCandidateSource() {
	}

	/**
	 * @param substitutions the reassigner's result for this validation; null
	 *                      gives no peak candidates
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
		if (substitutions != null) {
			for (SubstitutionReassigner.Proposal p : substitutions.proposals) {
				if (accepted.contains(p.psm) && agreesWithIdentification(p)) {
					result.candidates.add(candidate(p, true, undecidableToleranceDa));
					result.shapedCandidates++;
				}
			}
			for (SubstitutionReassigner.Proposal p : substitutions.rejectedPeakProposals) {
				if (accepted.contains(p.psm) && agreesWithIdentification(p)) {
					result.candidates.add(candidate(p, false, undecidableToleranceDa));
					result.rejectedPeakCandidates++;
				}
			}
		}
		addLocalised(result, accepted, undecidableToleranceDa);
		return result;
	}

	/**
	 * Whether a peak-level proposal is also the identification's call for its PSM. With the joint
	 * identification a PSM is a substitution candidate only as {@code OPEN_SUBSTITUTION_SITE}, at the
	 * residue {@link ShiftLocaliser} placed the shift on and for a residue its explanation names; a
	 * proposal made before scoring, from the search engine's placement, that says otherwise would look
	 * a peptide up that the identification does not claim. Any other identification (no group on the
	 * PSM) keeps every proposal, as before.
	 */
	static boolean agreesWithIdentification(SubstitutionReassigner.Proposal p) {
		OpenSearchPsm psm = p.psm;
		if (psm.getIdentificationGroup() == null) {
			return true;
		}
		ShiftLocaliser.Result loc = psm.getShiftLocalisation();
		if (!"OPEN_SUBSTITUTION_SITE".equals(psm.getIdentificationGroup()) || loc == null
				|| loc.state != ShiftLocaliser.State.LOCALISED_EXPLAINED || loc.site != p.position) {
			return false;
		}
		char to = p.substitution.to == 'I' ? 'L' : p.substitution.to;
		for (String label : loc.explanations) {
			if (label.startsWith(SUBSTITUTION_LABEL) && label.length() == SUBSTITUTION_LABEL.length() + 4) {
				char named = label.charAt(SUBSTITUTION_LABEL.length() + 3);
				if ((named == 'I' ? 'L' : named) == to) {
					return true;
				}
			}
		}
		return false;
	}

	/**
	 * One candidate per substitution in a validated PSM's explanation at its localised site. A peak
	 * proposal for the same PSM, site and residue gains the localisation instead of a second candidate.
	 */
	static void addLocalised(Result result, Set<OpenSearchPsm> accepted, double undecidableToleranceDa) {
		Map<String, VariantCandidate> byKey = new HashMap<>();
		for (VariantCandidate c : result.candidates) {
			byKey.put(key(c.file, c.scan, c.basePeptide, c.position, c.toAa), c);
		}
		List<OpenSearchPsm> psms = new ArrayList<>(accepted);
		psms.sort(Comparator.comparing(OpenSearchPsm::getFileName).thenComparingInt(OpenSearchPsm::getScan)
				.thenComparing(OpenSearchPsm::getPeptide));
		for (OpenSearchPsm psm : psms) {
			ShiftLocaliser.Result loc = psm.getShiftLocalisation();
			// Localised groups only: an established shift's explanation (deamidation = N->D) is the
			// peak's modification, not a variant call.
			if (loc == null || loc.state != ShiftLocaliser.State.LOCALISED_EXPLAINED
					|| !("OPEN_LOCALISED".equals(psm.getIdentificationGroup())
							|| "OPEN_SUBSTITUTION_SITE".equals(psm.getIdentificationGroup()))) {
				continue;
			}
			String peptide = psm.getPeptide();
			int position = loc.site + 1;
			for (String label : loc.explanations) {
				// "subst X->Y"
				if (!label.startsWith(SUBSTITUTION_LABEL) || label.length() != SUBSTITUTION_LABEL.length() + 4) {
					continue;
				}
				char from = label.charAt(SUBSTITUTION_LABEL.length());
				char to = label.charAt(SUBSTITUTION_LABEL.length() + 3);
				String k = key(psm.getFileName(), psm.getScan(), peptide, position, to);
				VariantCandidate existing = byKey.get(k);
				if (existing != null) {
					existing.localisationSiteProbability = loc.siteProbability;
					existing.localisationSupport = loc.support;
					continue;
				}
				VariantCandidate c = fromPsm(psm);
				c.position = position;
				c.fromAa = from;
				c.toAa = to;
				c.variantPeptide = peptide.substring(0, loc.site) + to + peptide.substring(loc.site + 1);
				c.deltaMass = SubstitutionMassTable.residueMass(to) - SubstitutionMassTable.residueMass(from);
				c.isobaricModifications.addAll(IsobaricModificationTable.sameResidue(from, c.deltaMass, position == 1,
						undecidableToleranceDa));
				c.ms2Decidable = c.isobaricModifications.isEmpty();
				c.frontEndEvidence = LOCALISED_CALL;
				c.localisationSiteProbability = loc.siteProbability;
				c.localisationSupport = loc.support;
				result.candidates.add(c);
				byKey.put(k, c);
				result.localisedCandidates++;
			}
		}
	}

	private static String key(String file, int scan, String peptide, int position, char to) {
		return file + '\t' + scan + '\t' + peptide + '\t' + position + '\t' + to;
	}

	static VariantCandidate candidate(SubstitutionReassigner.Proposal p, boolean peakShaped,
			double undecidableToleranceDa) {
		OpenSearchPsm psm = p.psm;
		VariantCandidate c = fromPsm(psm);
		c.position = p.position + 1;
		c.fromAa = p.substitution.from;
		c.toAa = p.substitution.to;
		c.variantPeptide = p.candidate;
		c.deltaMass = p.substitution.massDelta;
		c.peakSubstitutionShaped = peakShaped;
		c.isobaricModifications.addAll(IsobaricModificationTable.sameResidue(c.fromAa, c.deltaMass, c.position == 1,
				undecidableToleranceDa));
		c.ms2Decidable = c.isobaricModifications.isEmpty();
		c.frontEndEvidence = peakShaped ? String.valueOf(psm.getSubstitutionVerdict()) : "peak_not_substitution_shaped";
		c.reassignerInDb = p.inDatabase;
		return c;
	}

	/** The PSM-level fields every open-search candidate shares. */
	private static VariantCandidate fromPsm(OpenSearchPsm psm) {
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
		c.peakId = psm.getPeakId();
		c.localizedSequence = psm.getLocalizedSequence() == null ? "" : psm.getLocalizedSequence();
		return c;
	}

	static final double PROTON = 1.007276466812;
}
