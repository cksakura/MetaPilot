package uk.ac.quadram.metapilot.core.variant;

import java.util.Set;
import java.util.TreeSet;

/**
 * One amino-acid substitution candidate: one PSM x one substitution site.
 * <p>
 * The single record both front ends produce and everything downstream consumes
 * (lookup, control, classes, output files). A peptide with several
 * substitutions yields one record per site; all of them carry the same
 * {@link #variantPeptide}, with every substitution applied.
 */
public final class VariantCandidate {

	public enum Engine {
		PFIND, OPEN_SEARCH
	}

	/** Per-candidate class, identical vocabulary for both engines. */
	public enum PangenomeClass {
		/** Same-residue modification of equal mass; the pangenome columns are tie-break evidence only. */
		UNDECIDABLE_BY_MS2,
		/**
		 * Variant peptide is in a database this run searched: any sample's, not
		 * necessarily the one its own spectrum was searched against.
		 */
		SUBSTITUTION_IN_SEARCHED_DB,
		/** Not in the searched databases, found in UHGP-100. */
		SUBSTITUTION_IN_PANGENOME,
		/** In neither. */
		SUBSTITUTION_NOVEL
	}

	/**
	 * Outcome of the open-search MS2 form ranking ({@link VariantFormRanker}),
	 * among the forms at the candidate's own site.
	 */
	public enum Ms2Verdict {
		/** pFind: its call is trusted as-is, no ranking. */
		NOT_APPLICABLE,
		/** Open search, but no spectrum, no prediction or ranking switched off. */
		NOT_RANKED,
		/** This candidate's substitution is the best same-site form by more than the margin. */
		SUBSTITUTION,
		/** Another same-mass substitution of the same residue. */
		OTHER_SUBSTITUTION,
		/** The delta mass as a modification of the same residue. */
		MODIFICATION,
		AMBIGUOUS
	}

	public Engine engine;
	public String file;
	public int scan;
	/**
	 * Third part of the join key back to the workflow's PSM table: the hit rank
	 * for open search ({@code validated_psms.tsv}, {@code substitution_forms.tsv}),
	 * the precursor charge for pFind ({@code combined/psms.tsv}).
	 */
	public int psmId;
	/** Sequence as identified; I and L kept as reported. */
	public String basePeptide;
	/** 1-based residue index in {@link #basePeptide}. */
	public int position;
	public char fromAa;
	public char toAa;
	public String variantPeptide;
	/** {@code mass(to) - mass(from)}, underivatised, with the engine's cysteine convention. */
	public double deltaMass;
	/** Open search: the delta-mass peak; pFind: -1. */
	public int peakId = -1;
	/**
	 * Open search: whether {@code SubstitutionReassigner} judged the peak
	 * substitution-shaped. Null for pFind, which has no peaks.
	 */
	public Boolean peakSubstitutionShaped;
	public boolean ms2Decidable = true;
	/** Same-residue modifications of equal mass, or pFind's own isobaric alias. */
	public final Set<String> isobaricModifications = new TreeSet<>();
	/**
	 * What the front end called: pFind's raw modification token; open search,
	 * the reassigner's verdict, or {@code peak_not_substitution_shaped}.
	 */
	public String frontEndEvidence = "";
	/** Open search: the engine's localized sequence; empty for pFind. */
	public String localizedSequence = "";
	/** Open search: whether the reassigner found the variant in its database; null for pFind. */
	public Boolean reassignerInDb;
	/** Precursor charge; 0 when unknown. */
	public int charge;
	/** Observed precursor m/z; NaN when unknown. */
	public double precursorMz = Double.NaN;
	/** The PSM's q-value as the workflow reported it; NaN when unknown. */
	public double qValue = Double.NaN;
	/** Searched-database proteins the PSM was assigned to, as the engine reported them. */
	public final Set<String> assignedProteins = new TreeSet<>();
	/**
	 * Where {@link #basePeptide} sits in each assigned protein it was found in
	 * among the searched databases; filled by the lookup.
	 */
	final java.util.Map<String, PangenomeVariantLookup.Location> proteinLocations = new java.util.TreeMap<>();
	/**
	 * Candidates from the same PSM, this one included: pFind's sites of one
	 * multi-substitution peptide, or the open-search reassigner's alternative
	 * readings of one delta mass.
	 */
	public int candidatesForPsm = 1;

	// --- filled by the MS2 form ranking (open search) ---------------------
	// Secondary evidence: it does not set pangenomeClass.
	public Ms2Verdict ms2Verdict = Ms2Verdict.NOT_APPLICABLE;
	/** Spectral-entropy similarity of this candidate's substituted peptide; NaN when not ranked. */
	public double ms2ScoreSubstitution = Double.NaN;
	/** Best other same-mass substitution of the same residue; NaN when there is none or not ranked. */
	public double ms2ScoreOtherSubstitution = Double.NaN;
	/** The delta mass as a modification of the substituted residue. */
	public double ms2ScoreModification = Double.NaN;
	/** Best form changing any other residue, substitution or modification; NaN when none was scored. */
	public double ms2ScoreElsewhere = Double.NaN;
	/**
	 * Whether a form at another residue beats every same-site form by more than
	 * the margin: the delta mass seems to sit elsewhere. Null when not ranked.
	 */
	public Boolean ms2LocalisesElsewhere;
	/** Best-scoring form at any site. */
	public String ms2BestForm = "";

	// --- filled by the lookup -------------------------------------------
	/** Null when no searched database could be read. */
	public Boolean inSearchedDb;
	public boolean inPangenomeAll;
	/** Null when DETECTED scope is unavailable or not reported. */
	public Boolean inPangenomeDetected;
	/**
	 * Whether UHGP-100 has the variant in a species the PSM's own proteins
	 * belong to: natural variation of that species, rather than a sequence some
	 * unrelated gut microbe happens to carry. Null when the assigned proteins'
	 * species or the hits' species are unknown.
	 */
	public Boolean inSameSpecies;
	/** Species representatives of {@link #assignedProteins}. */
	public final Set<String> assignedSpecies = new TreeSet<>();
	/** -1 when species could not be resolved. */
	public int nSpecies = -1;
	public String exampleProtein = "";

	// --- filled from the UHGP-100 cluster membership; -1 / "" when unavailable ---
	/** Genomes, any species, with a protein in a UHGP-100 cluster that has the variant. */
	public int nGenomes = -1;
	/** Of {@link #nGenomes}, those in {@link #assignedSpecies}. */
	public int nGenomesSameSpecies = -1;
	/** Every catalogued genome of {@link #assignedSpecies}, carrying the variant or not. */
	public int nGenomesInSpecies = -1;
	/** Of {@link #nGenomesSameSpecies}, isolate (cultured) genomes rather than MAGs. */
	public int nIsolatesSameSpecies = -1;
	/** A same-species genome with the variant: the lowest-numbered isolate, else the lowest-numbered genome. */
	public String exampleGenomeSameSpecies = "";

	// --- eggNOG function of assignedProteins, union over them; null (NA) when no table covers them ---
	public String assignedKeggKo;
	public String assignedCog;
	public String assignedCogCategory;
	public String assignedPreferredName;
	public String assignedDescription;

	public PangenomeClass pangenomeClass;

	/** Label used for grouping: I and L are one residue to a mass spectrometer. */
	public String substitutionLabel() {
		return residueLabel(fromAa) + "->" + residueLabel(toAa);
	}

	static String residueLabel(char aa) {
		return aa == 'I' || aa == 'L' ? "I/L" : String.valueOf(aa);
	}

	/** Residues where {@link #variantPeptide} differs from {@link #basePeptide}, I and L counted as one. */
	public int substitutionsInVariantPeptide() {
		if (basePeptide == null || variantPeptide == null || basePeptide.length() != variantPeptide.length()) {
			return -1;
		}
		String a = PeptideAutomaton.normalise(basePeptide);
		String b = PeptideAutomaton.normalise(variantPeptide);
		int n = 0;
		for (int i = 0; i < a.length(); i++) {
			if (a.charAt(i) != b.charAt(i)) {
				n++;
			}
		}
		return n;
	}

	/** 1-based protein index of the substituted residue, or -1 when the peptide was not located there. */
	public int proteinPosition(String protein) {
		PangenomeVariantLookup.Location at = proteinLocations.get(protein);
		return at == null ? -1 : at.start + position - 1;
	}

	/**
	 * Undecidable takes precedence over every database answer: a same-residue
	 * modification of equal mass explains the spectrum as well as a database hit
	 * does, so a hit cannot promote the candidate to a substitution.
	 */
	void classify() {
		if (!ms2Decidable) {
			pangenomeClass = PangenomeClass.UNDECIDABLE_BY_MS2;
		} else if (Boolean.TRUE.equals(inSearchedDb)) {
			pangenomeClass = PangenomeClass.SUBSTITUTION_IN_SEARCHED_DB;
		} else if (inPangenomeAll) {
			pangenomeClass = PangenomeClass.SUBSTITUTION_IN_PANGENOME;
		} else {
			pangenomeClass = PangenomeClass.SUBSTITUTION_NOVEL;
		}
	}
}
