package uk.ac.quadram.metapilot.core.opensearch;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.Set;
import java.util.function.Consumer;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

/**
 * Decides, for every PSM sitting on a delta mass that a single amino-acid
 * substitution could explain, whether it is a modification or a misassigned
 * sequence - and proposes the sequence it should have been.
 *
 * <h2>The problem this exists for</h2>
 * Measured on 12mix (2026-09-14). The two accepted peaks at +/-14.0157 Da (CH2)
 * passed the per-peak local-FDR walk at 0.94% and 0.96% decoy, and their
 * <em>true</em> error rate measured against the twelve spiked genomes was 18.35%
 * and 17.85% - a nineteen-fold under-estimate. +28.0312 (two CH2) was 14.10%.
 * Every other accepted peak sat at or below 4.56%.
 * <p>
 * Nothing else in this workflow can see it:
 * <ul>
 * <li>Target-decoy is blind by construction - the wrong peptide is itself a
 * target sequence, so no reversed or shuffled decoy ever competes with it.</li>
 * <li>{@link DeltaMassEvidenceGate} passes these peaks, because a substituted
 * residue genuinely shifts part of the fragment series by the substitution mass.
 * The evidence is real; the interpretation is wrong.</li>
 * <li>Peptide corroboration passes them too - the true peptide really is in the
 * zero family, which is how the substitution was detected in the first
 * place.</li>
 * </ul>
 *
 * <h2>What the data says the method must handle</h2>
 * Generating every CH2 substitution of the 2,369 PSMs on those two peaks and
 * looking each up in the full tryptic space of the searched databases
 * (6,532,401 distinct peptides) gave:
 * <pre>
 *   substituted peptide in DB and already identified unmodified    601   25.4%
 *   substituted peptide in DB, not identified                       98    4.1%
 *   substituted peptide NOT in the searched database             1,670   70.5%
 * </pre>
 * Entrapment of the assigned peptide was 20.41% / 18.37% / 18.44% - uniform, so
 * one rule covers all three. <b>The dominant case is a search-space gap, not a
 * scoring error.</b> Seven times in ten the true sequence was never a candidate,
 * because genome selection picked one strain and the peptide came from a
 * relative that differs by one residue. That is specific to metaproteomics, and
 * it is why this class proposes candidates for a re-search rather than simply
 * preferring an already-searched alternative.
 *
 * <h2>Localization is what resolves the mass collisions</h2>
 * Several substitution masses coincide with real modifications - N-&gt;D with
 * deamidation, A-&gt;S and F-&gt;Y with oxidation, A-&gt;Q and G-&gt;N with
 * carbamidomethyl. On the same run those peaks carried 22,186 / 4,790 / 1,913
 * PSMs at 2.80% / 1.68% / 1.17% entrapment: overwhelmingly the modification. A
 * mass match alone would condemn them. So a candidate is only generated when the
 * delta <em>localizes to a residue the substitution could act on</em>: if
 * +15.9949 sits on an M it is oxidation, and A-&gt;S is never proposed. This
 * also collapses the candidate count per PSM from roughly nine to usually one.
 *
 * <h2>What this class does and does not do</h2>
 * It classifies, flags and proposes. It does <b>not</b> rewrite a PSM's peptide,
 * because choosing between the shifted assignment and the substituted one
 * requires scoring the spectrum against the candidate, and this workflow has no
 * scorer - MSFragger is the scorer. {@link #writeCandidateFasta} emits the
 * proposals as a FASTA for a bounded second search; that pass is what can settle
 * the 70.5% whose sequence is absent from the database, which no lookup can.
 * Until that pass runs, a substitution-shaped peak is reported as such and kept
 * out of the modification verdict rather than counted as a PTM.
 */
public final class SubstitutionReassigner {

	private static final Logger LOGGER = LogManager.getLogger(SubstitutionReassigner.class);

	/** How a PSM on a substitution-shaped peak was classified. */
	public enum Verdict {
		/** Delta does not localize to a residue any candidate substitution fits. */
		NOT_SUBSTITUTION,
		/** Candidate exists in the searched database; a re-search can confirm it. */
		SUBSTITUTION_IN_DB,
		/** Candidate is absent from the searched database; needs the augmented pass. */
		SUBSTITUTION_NOT_IN_DB
	}

	public static final class Config {
		/**
		 * Tolerance for matching a peak centroid to a substitution mass. Fitted
		 * peaks have sigma near 2 mDa, so 5 mDa admits the peak and little else.
		 */
		public double massToleranceDa = 0.005;
		/**
		 * A peak is declared substitution-shaped when at least this share of its
		 * members yield a localization-consistent candidate. Oxidation sits far
		 * below it because its delta localizes to residues the colliding
		 * substitutions cannot act on; deamidation is handled by
		 * {@link #modificationEquivalent}.
		 */
		public double minSubstitutionFraction = 0.50;
		/** Below this many members the fraction is too noisy to act on. */
		public int minPsms = 50;
		/**
		 * Minimum share of a peak's candidates that must exist in the searched
		 * database before it is called a substitution rather than a modification.
		 * <p>
		 * This is the criterion that separates the two, and localization is not.
		 * Measured on 12mix (2026-09-14), on the first build that ran without it:
		 * localization alone withdrew 22 peaks, of which the three largest were
		 * deamidation (22,019 PSMs, 2.60% entrapment), +40.0313 (6,454, 1.63%) and
		 * carbamidomethyl (2,011, 1.16%) - real modifications, wrongly condemned.
		 * The reason localization cannot help there is structural: deamidation acts
		 * on N and Q, and so do N-&gt;D and Q-&gt;E; +40.0313 acts on G and so does
		 * G-&gt;P. The substitution acts on exactly the residue the modification
		 * does, so the position agrees with both.
		 * <p>
		 * Database membership does separate them, because the two arise
		 * differently. A real substitution exists because the true peptide is in a
		 * <em>related genome</em>, and in a metaproteome whose database is a set of
		 * selected strains that partner is often present. A real modification's
		 * "substituted" candidate is an arbitrary sequence that essentially never
		 * is. Sorted by that share, the same 22 peaks split at a gap:
		 * <pre>
		 *   0.00% +40.0313 | 0.48% +57.0214 | 1.76% +0.9841 | 3.48% -14.9698 | 4.66% +58.0056
		 *   -- gap --
		 *   6.54% -58.0052 | 9.09% +28.0312 | 16.82% +14.0156 | 19.08% -14.0157 | 21.50% +27.0101
		 *
		 *   below 5%:  5 peaks, 30,828 PSMs, weighted entrapment  2.38%
		 *   at/above:  11 peaks, 4,618 PSMs, weighted entrapment 15.33%
		 * </pre>
		 * A 6.4-fold separation. One peak is misfiled at this threshold, +58.0056
		 * at 4.66% in-DB and 12.32% entrapment, costing 224 PSMs - the price of a
		 * rule that keeps 30,604 correct modification PSMs.
		 * <p>
		 * Set to 0 to disable the test, which restores the localization-only
		 * behaviour and its false positives.
		 */
		public double minInDatabaseFraction = 0.05;
		/**
		 * Substitutions that yield the same molecule as a common modification of
		 * the same residue: N-&gt;D and Q-&gt;E are deamidated N and Q. Mass,
		 * localization and fragments are identical, so MS/MS cannot separate them,
		 * and the modification is always the reported call. The in-database test
		 * does not hold them apart either: on PXD050781 (MSFragger open,
		 * 2026-10-12) 10.2% of the deamidation peak's N-&gt;D/Q-&gt;E candidates
		 * existed in a database holding E. coli and several Blautia relatives,
		 * above the 5% threshold, and 19,485 deamidated PSMs were withdrawn.
		 */
		public Set<String> modificationEquivalent = Set.of("N->D", "Q->E");
		/** Off restores the previous behaviour exactly. */
		public boolean enabled = true;
	}

	/** One proposal: this PSM's spectrum probably came from {@code candidate}. */
	public static final class Proposal {
		public final OpenSearchPsm psm;
		public final String candidate;
		public final SubstitutionMassTable.Substitution substitution;
		/** 0-based position in the assigned peptide. */
		public final int position;
		public final boolean inDatabase;

		Proposal(OpenSearchPsm psm, String candidate, SubstitutionMassTable.Substitution substitution, int position,
				boolean inDatabase) {
			this.psm = psm;
			this.candidate = candidate;
			this.substitution = substitution;
			this.position = position;
			this.inDatabase = inDatabase;
		}

		/** Human-readable change, 1-based as sequence positions are normally cited. */
		public String describe() {
			return String.format(Locale.ROOT, "%c%d%c", substitution.from, position + 1, substitution.to);
		}
	}

	public static final class PeakSummary {
		public final int peakId;
		public final double centroid;
		public final int psms;
		public final int withCandidate;
		public final int candidateInDb;
		public final boolean substitutionShaped;
		public final List<String> substitutions;

		PeakSummary(int peakId, double centroid, int psms, int withCandidate, int candidateInDb,
				boolean substitutionShaped, List<String> substitutions) {
			this.peakId = peakId;
			this.centroid = centroid;
			this.psms = psms;
			this.withCandidate = withCandidate;
			this.candidateInDb = candidateInDb;
			this.substitutionShaped = substitutionShaped;
			this.substitutions = substitutions;
		}

		public double substitutionFraction() {
			return psms == 0 ? 0.0 : (double) withCandidate / psms;
		}

		/**
		 * Share of this peak's candidates that exist in the searched database -
		 * the criterion that separates a substitution from a modification whose
		 * mass happens to match one. See {@link Config#minInDatabaseFraction}.
		 */
		public double inDatabaseFraction() {
			return withCandidate == 0 ? 0.0 : (double) candidateInDb / withCandidate;
		}
	}

	public static final class Result {
		public final List<Proposal> proposals;
		public final List<PeakSummary> peaks;
		/** Peaks to keep out of the modification verdict. */
		public final Set<Integer> substitutionPeakIds;
		/**
		 * Proposals on peaks whose centroid matches a substitution mass but which
		 * were not judged substitution-shaped. Report-only: no verdict, candidate
		 * or reassignment is set on these PSMs, and nothing in this class or the
		 * workflow acts on them. They exist so a downstream report can compare the
		 * modification-like peaks (+0.9841, +40.0313, +57.0214, ...) with the
		 * substitution-shaped ones on the same footing.
		 */
		public final List<Proposal> rejectedPeakProposals;

		Result(List<Proposal> proposals, List<PeakSummary> peaks, Set<Integer> substitutionPeakIds,
				List<Proposal> rejectedPeakProposals) {
			this.proposals = proposals;
			this.peaks = peaks;
			this.substitutionPeakIds = substitutionPeakIds;
			this.rejectedPeakProposals = rejectedPeakProposals;
		}

		public int inDatabaseCount() {
			int n = 0;
			for (Proposal p : proposals) {
				if (p.inDatabase) {
					n++;
				}
			}
			return n;
		}
	}

	private SubstitutionReassigner() {
	}

	/**
	 * Classifies every PSM on a substitution-shaped peak and marks the ones that
	 * are misassigned sequences rather than modifications.
	 *
	 * @param psms  every candidate, mutated in place where a verdict is reached
	 * @param index peptide space of the databases actually searched; null skips
	 *              the database test and reports every candidate as absent
	 */
	public static Result analyse(List<OpenSearchPsm> psms, TrypticPeptideIndex index, Config cfg,
			Consumer<String> progress) {
		if (cfg == null) {
			cfg = new Config();
		}
		if (!cfg.enabled) {
			return new Result(new ArrayList<>(), new ArrayList<>(), new HashSet<>(), new ArrayList<>());
		}

		// Group by peak so the decision is made per mass hypothesis, not per PSM:
		// a single PSM whose delta happens to fit a substitution proves nothing.
		Map<Integer, List<OpenSearchPsm>> byPeak = new HashMap<>();
		Map<Integer, Double> centroid = new HashMap<>();
		for (OpenSearchPsm psm : psms) {
			int pid = psm.getPeakId();
			if (pid < 0) {
				continue;
			}
			double c = psm.getPeakCentroid();
			if (Double.isNaN(c) || Math.abs(c) < cfg.massToleranceDa) {
				continue;
			}
			centroid.putIfAbsent(pid, c);
			byPeak.computeIfAbsent(pid, k -> new ArrayList<>()).add(psm);
		}

		List<Proposal> proposals = new ArrayList<>();
		List<Proposal> rejectedPeakProposals = new ArrayList<>();
		List<PeakSummary> summaries = new ArrayList<>();
		Set<Integer> substitutionPeaks = new HashSet<>();

		for (Map.Entry<Integer, List<OpenSearchPsm>> e : byPeak.entrySet()) {
			int pid = e.getKey();
			double c = centroid.get(pid);
			List<SubstitutionMassTable.Substitution> subs = new ArrayList<>(
					SubstitutionMassTable.matching(c, cfg.massToleranceDa));
			if (cfg.modificationEquivalent != null) {
				Set<String> equivalent = cfg.modificationEquivalent;
				subs.removeIf(s -> equivalent.contains(s.from + "->" + s.to));
			}
			if (subs.isEmpty()) {
				continue;
			}
			List<Proposal> peakProposals = new ArrayList<>();
			int inDb = 0;
			for (OpenSearchPsm psm : e.getValue()) {
				Proposal p = propose(psm, subs, index);
				if (p != null) {
					peakProposals.add(p);
					if (p.inDatabase) {
						inDb++;
					}
				}
			}
			int n = e.getValue().size();
			double frac = n == 0 ? 0.0 : (double) peakProposals.size() / n;
			double inDbShare = peakProposals.isEmpty() ? 0.0 : (double) inDb / peakProposals.size();
			// Both criteria, not either: the delta must be explainable by a
			// substitution the localization agrees with, AND those substituted
			// sequences must actually turn up in the database often enough to be
			// coming from a relative rather than from arithmetic.
			boolean shaped = n >= cfg.minPsms && frac >= cfg.minSubstitutionFraction
					&& inDbShare >= cfg.minInDatabaseFraction;
			Set<String> names = new LinkedHashSet<>();
			for (SubstitutionMassTable.Substitution s : subs) {
				names.add(s.from + "->" + s.to);
			}
			summaries.add(new PeakSummary(pid, c, n, peakProposals.size(), inDb, shaped, new ArrayList<>(names)));
			if (shaped) {
				substitutionPeaks.add(pid);
				proposals.addAll(peakProposals);
				for (Proposal p : peakProposals) {
					p.psm.setSubstitutionVerdict(
							(p.inDatabase ? Verdict.SUBSTITUTION_IN_DB : Verdict.SUBSTITUTION_NOT_IN_DB).name());
					p.psm.setSubstitutionCandidate(p.candidate);
					p.psm.setDeltaMassReassignment("SUBSTITUTION:" + p.describe());
				}
			} else {
				rejectedPeakProposals.addAll(peakProposals);
			}
		}
		summaries.sort((a, b) -> Double.compare(b.substitutionFraction(), a.substitutionFraction()));
		Result result = new Result(proposals, summaries, substitutionPeaks, rejectedPeakProposals);

		if (progress != null) {
			progress.accept(String.format(Locale.ROOT,
					"Substitution analysis: %d of %d non-zero peaks are substitution-shaped (>=%.0f%% of members "
							+ "yield a localization-consistent single-residue candidate); %,d PSMs proposed, "
							+ "%,d of them with the candidate already in the searched database",
					substitutionPeaks.size(), byPeak.size(), 100 * cfg.minSubstitutionFraction, proposals.size(),
					result.inDatabaseCount()));
			int shown = 0;
			for (PeakSummary s : summaries) {
				if (!s.substitutionShaped || shown++ >= 8) {
					continue;
				}
				progress.accept(String.format(Locale.ROOT,
						"    peak %d at %+.4f Da (%s): %.1f%% substitution-consistent of %d, %.1f%% of candidates in DB",
						s.peakId, s.centroid, String.join(", ", s.substitutions), 100 * s.substitutionFraction(),
						s.psms, 100 * s.inDatabaseFraction()));
			}
		}
		LOGGER.info("Substitution analysis: {} substitution-shaped peaks, {} proposals, {} with the candidate in the "
				+ "searched database", substitutionPeaks.size(), proposals.size(), result.inDatabaseCount());
		return result;
	}

	/**
	 * The one candidate this PSM's localization supports, or null.
	 * <p>
	 * Only positions the delta actually localizes to are considered. MSFragger
	 * marks them by lower-casing the residue in
	 * {@link OpenSearchPsm#getLocalizedSequence()}; when it localized nothing -
	 * ambiguous placement, or no localization payload at all - every position is
	 * allowed, which is the conservative reading: it lets a candidate be proposed
	 * for auditing rather than silently declaring the PSM a modification.
	 * <p>
	 * When more than one position fits, the one already present in the database
	 * wins; failing that, the first. A PSM with several equally good candidates is
	 * exactly the case the re-search is for.
	 */
	static Proposal propose(OpenSearchPsm psm, List<SubstitutionMassTable.Substitution> subs,
			TrypticPeptideIndex index) {
		String pep = psm.getPeptide();
		if (pep == null || pep.isEmpty()) {
			return null;
		}
		boolean[] localized = localizedPositions(psm, pep.length());
		Proposal first = null;
		for (int i = 0; i < pep.length(); i++) {
			if (!localized[i]) {
				continue;
			}
			char at = pep.charAt(i);
			for (SubstitutionMassTable.Substitution s : subs) {
				if (s.from != at) {
					continue;
				}
				String candidate = pep.substring(0, i) + s.to + pep.substring(i + 1);
				boolean inDb = index != null && index.contains(candidate);
				Proposal p = new Proposal(psm, candidate, s, i, inDb);
				if (inDb) {
					return p;
				}
				if (first == null) {
					first = p;
				}
			}
		}
		return first;
	}

	/** Positions the delta was localized to, or all positions when it was not. */
	private static boolean[] localizedPositions(OpenSearchPsm psm, int length) {
		boolean[] out = new boolean[length];
		String loc = psm.getLocalizedSequence();
		int marked = 0;
		if (loc != null && loc.length() == length) {
			for (int i = 0; i < length; i++) {
				if (Character.isLowerCase(loc.charAt(i))) {
					out[i] = true;
					marked++;
				}
			}
		}
		if (marked == 0) {
			java.util.Arrays.fill(out, true);
		}
		return out;
	}

	/**
	 * Writes the proposed sequences as a FASTA for a bounded second MSFragger
	 * pass.
	 * <p>
	 * This is the interface to the only step that can settle the 70.5% of cases
	 * whose true sequence is not in any searched database: a lookup cannot confirm
	 * a peptide that was never there, but a search against these generated
	 * sequences can. One entry per distinct candidate, named after the assigned
	 * peptide and the substitution so a hit is traceable back to the PSM that
	 * proposed it.
	 *
	 * @return how many distinct candidate sequences were written
	 */
	public static int writeCandidateFasta(Result result, File out) throws IOException {
		out.getParentFile().mkdirs();
		Map<String, String> byCandidate = new java.util.LinkedHashMap<>();
		for (Proposal p : result.proposals) {
			byCandidate.putIfAbsent(p.candidate,
					"sub_" + p.psm.getPeptide() + "_" + p.describe() + (p.inDatabase ? "_inDB" : ""));
		}
		try (PrintWriter w = new PrintWriter(out, "UTF-8")) {
			w.println("# Substitution candidates proposed by MetaPilot's open-search validation.");
			w.println("# Each entry is one assigned peptide with a single residue changed, where the delta mass");
			w.println("# and its localization are both consistent with that change. Search these against the");
			w.println("# affected spectra to decide whether the PSM is a modification or a misassigned sequence.");
			for (Map.Entry<String, String> e : byCandidate.entrySet()) {
				w.println(">" + e.getValue());
				w.println(e.getKey());
			}
		}
		return byCandidate.size();
	}

	/** Per-peak audit of the substitution analysis. */
	public static void writePeakReport(Result result, File out) throws IOException {
		out.getParentFile().mkdirs();
		try (PrintWriter w = new PrintWriter(out, "UTF-8")) {
			w.println("# A peak is substitution-shaped when most of its PSMs yield a single-residue candidate whose");
			w.println("# position agrees with where MSFragger localized the delta. Mass match alone is not enough:");
			w.println("# A->S and F->Y collide with oxidation, A->Q and G->N with carbamidomethyl; those peaks score");
			w.println("# low here because their delta localizes elsewhere. N->D and Q->E are not tested: they are");
			w.println("# the same molecule as deamidated N and Q, so the deamidation peak is always a modification.");
			w.println("peak_id\tcentroid_da\tpsms\twith_candidate\tsubstitution_fraction\tcandidate_in_db\t"
					+ "in_db_fraction\tsubstitution_shaped\tsubstitutions");
			for (PeakSummary s : result.peaks) {
				w.println(s.peakId + "\t" + String.format(Locale.ROOT, "%.5f", s.centroid) + "\t" + s.psms + "\t"
						+ s.withCandidate + "\t" + String.format(Locale.ROOT, "%.4f", s.substitutionFraction()) + "\t"
						+ s.candidateInDb + "\t" + String.format(Locale.ROOT, "%.4f", s.inDatabaseFraction()) + "\t"
						+ s.substitutionShaped + "\t" + String.join(",", s.substitutions));
			}
		}
	}
}
