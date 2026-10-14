package uk.ac.quadram.metapilot.dbSearch.sage;

import java.util.Locale;

/**
 * The two decoy populations a Sage MAG run produces, and which file each one
 * lives in. They are deliberately different, and conflating them breaks FDR in
 * opposite directions at the two levels - so the distinction is named here once
 * rather than left as two unexplained string literals in two readers.
 *
 * <h2>results.sage.tsv - Sage's own decoys</h2>
 * Sage generates these itself ({@code generate_decoys: true}) and tags them with
 * its configured {@code decoy_tag}. They are what Sage's own q-values are
 * estimated against, and the identification output retains them, so PSM-level
 * target/decoy comes straight from this tag. {@link SagePsmReader} uses it.
 *
 * <h2>lfq.tsv - the shuffled entries that survive quantification</h2>
 * Sage's LFQ stage quantifies only targets ({@code peptide_q <= 0.01 &&
 * label == 1}), exactly as Philosopher strips {@code rev_} from FragPipe's
 * quantified list. Genome-level FDR is computed after quantification, so it
 * would otherwise have no decoy population at all. The workflow therefore writes
 * its own {@code shuf_} sequences into the combined FASTA as ordinary entries:
 * Sage sees them as targets, quantifies them, and they reach lfq.tsv, where the
 * genome-level FDR uses them. {@link SageResultReader} uses this tag.
 *
 * <p>
 * So a {@code shuf_} hit is <em>not</em> a decoy to Sage, and a Sage-generated
 * decoy never reaches the quantified list. Each reader asks the question its own
 * file can answer.
 */
public final class SageDecoyTags {

	/** Prefixes Sage itself puts on the decoys it generates. */
	private static final String[] ENGINE_TAGS = { "rev_", "decoy_" };

	/** Prefix this application writes on the entries that carry genome-level FDR. */
	public static final String SHUFFLED_TAG = "shuf_";

	private SageDecoyTags() {
	}

	/**
	 * True for a decoy Sage generated - the target/decoy label of the
	 * identification output. Case-insensitive, because the configured tag and the
	 * tag the readers were written against have differed in case.
	 */
	public static boolean isEngineDecoy(String accession) {
		if (accession == null || accession.isEmpty()) {
			return false;
		}
		String lower = accession.toLowerCase(Locale.ROOT);
		for (String tag : ENGINE_TAGS) {
			if (lower.startsWith(tag)) {
				return true;
			}
		}
		return false;
	}

	/** True when any accession in {@code accessions} is a Sage-generated decoy. */
	public static boolean anyEngineDecoy(String[] accessions) {
		return anyMatch(accessions, true);
	}

	/**
	 * True for one of the shuffled entries this application wrote into the
	 * combined FASTA - the decoy population for genome-level FDR, which Sage
	 * scores as a target.
	 */
	public static boolean isShuffledDecoy(String accession) {
		return accession != null && accession.toLowerCase(Locale.ROOT).startsWith(SHUFFLED_TAG);
	}

	/** True when any accession in {@code accessions} is a shuffled entry. */
	public static boolean anyShuffledDecoy(String[] accessions) {
		return anyMatch(accessions, false);
	}

	private static boolean anyMatch(String[] accessions, boolean engine) {
		if (accessions == null) {
			return false;
		}
		for (String accession : accessions) {
			String trimmed = accession == null ? null : accession.trim();
			if (engine ? isEngineDecoy(trimmed) : isShuffledDecoy(trimmed)) {
				return true;
			}
		}
		return false;
	}
}
