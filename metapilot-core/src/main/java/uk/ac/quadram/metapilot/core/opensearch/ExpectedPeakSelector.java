package uk.ac.quadram.metapilot.core.opensearch;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.IdentityHashMap;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.Set;
import java.util.function.Consumer;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

/**
 * Stage 1, second pass: admits PSMs carrying an <em>expected</em> mass shift into
 * the genome-mapping evidence set, alongside the zero family that
 * {@link UnmodifiedPsmSelector} produces.
 *
 * <h2>Why this exists</h2>
 * Stage 1 restricted itself to the zero peak on the reasoning that "the rest of
 * the delta-mass axis is of no interest" to genome selection. That reasoning
 * fails in one specific way. A peptide carrying a mass shift maps to exactly the
 * same genome as its unmodified form, so the modification state is irrelevant to
 * genome selection — while the closed-search MAG workflows, which declare
 * oxidation and N-terminal acetyl as variable modifications, match those spectra
 * with {@code massdiff ~ 0} and keep them. Measured on the 12mix mock community
 * against its 12 known genomes, Met oxidation alone accounted for a fifth of the
 * PSMs this workflow was discarding, and the discarded set was 96% consistent
 * with the ground truth — as good as the zero family itself.
 *
 * <h2>Why it is not peak discovery</h2>
 * {@link UnmodifiedPsmSelector} declines to run the classifier at stage 1
 * because a well-trained classifier wants the refined database that stage 1
 * exists to produce. The same objection would apply to discovering mass shifts
 * here, so this class does not discover: it tests a fixed, auditable allow-list
 * of shifts that are expected a priori, and every one of them must earn its place
 * in <em>this run's</em> data. Discovery remains stage 2's job.
 *
 * <h2>The four criteria</h2>
 * A shift is admitted only if all of these hold, evaluated per run:
 * <ol>
 * <li><b>On the allow-list</b> ({@link #DEFAULT_SHIFTS}), within
 * {@link Config#toleranceDa} of the run's calibrated position. Assignment is
 * exclusive: the zero family claims its own members (including its isotope
 * satellites) first, and each remaining PSM goes to its nearest single shift.</li>
 * <li><b>The peak is present</b>: at least {@link Config#minTargetsInWindow}
 * target PSMs in the window. Below that the window cannot be tested and cannot
 * matter.</li>
 * <li><b>Target enrichment</b>: raw decoys/targets in the window at or below
 * {@link Config#maxWindowFdr}, before any cut. A window sitting near 1:1 is
 * background, not a peak. Measured separation on the 12mix: admitted shifts
 * 8.6-25.8%, rejected ones 57.9-88.6%.</li>
 * <li><b>Peptide corroboration</b>: at least {@link Config#minCorroboration} of
 * the shift's accepted target peptides were <em>also</em> accepted unmodified in
 * the same run's zero family. This is the decisive criterion, and it is the one
 * that keeps the step inside stage 1's remit. A modified peptide whose
 * unmodified form is already confidently identified is corroborated evidence for
 * a genome already implied; a "modified peptide" never seen unmodified is an
 * unsupported new claim, which stage 1 has no business making. Measured on the
 * 12mix, every shift whose accepted set was &gt;=83% ground-truth consistent had
 * corroboration 59.8-94.5%, and every shift below 38% had 0-14.3% — the two
 * populations do not overlap.</li>
 * </ol>
 * Each admitted shift then gets its own q-value cut, computed <em>within</em>
 * that peak and never pooled, so a healthy zero peak cannot carry a sick
 * oxidation peak past the threshold.
 * <p>
 * Finally {@link Config#maxAdmittedFraction} caps the total admitted against the
 * zero family's own accepted target count, so a pathological run cannot end up
 * dominated by shifted evidence. On the 12mix the real figure was ~5%.
 *
 * <h2>Measured effect</h2>
 * Nine 12mix runs, marker stage: +2,877 accepted PSMs on 52,970 (+5.4%) at 94.6%
 * entrapment consistency, against the zero family's own 93.9%. The change is
 * purity-positive, not a trade.
 *
 * <h2>Known limitation</h2>
 * Isotope satellites of a shifted peak are not folded onto their root here, only
 * those of the zero family are (which {@link UnmodifiedPsmSelector} does before
 * this class runs). A +1 satellite of the oxidation peak is therefore left
 * unassigned rather than admitted. That is a yield question, not a correctness
 * one, and folding it would need the same consecutive-satellite and
 * amplitude-monotonicity tests stage 2 applies before it can be trusted.
 */
public final class ExpectedPeakSelector {

	private static final Logger LOGGER = LogManager.getLogger(ExpectedPeakSelector.class);

	/** One entry on the allow-list: a name for the audit trail and a mass in Da. */
	public static final class ExpectedShift {
		public final String name;
		public final double massDa;

		public ExpectedShift(String name, double massDa) {
			this.name = name;
			this.massDa = massDa;
		}

		@Override
		public String toString() {
			return String.format(Locale.ROOT, "%s (%+.4f Da)", name, massDa);
		}
	}

	/**
	 * The allow-list, as measured on the 12mix mock community against its 12 known
	 * genomes. Every entry here passed all four criteria in all nine runs; the
	 * candidates that failed are listed below so the list is auditable in both
	 * directions.
	 * <p>
	 * Note that these are <b>expected mass shifts</b>, not strictly biological
	 * modifications: Na-H and Ca-2H are metal adducts and carbamidomethyl here is
	 * over-alkylation beyond the fixed modification. For genome mapping that
	 * distinction does not matter — they are real peptides — and dropping them
	 * would discard two of the cleanest entries measured (96.8% and 95.1%
	 * ground-truth consistency).
	 * <p>
	 * <b>Deliberately excluded, with measured reasons</b> (12mix, nine runs,
	 * ground-truth consistency of the accepted set / peptide corroboration):
	 * methyl +14.0157 (37.9% / 12.1%), dimethyl +28.0313 (33.3% / 6.3%), Lys loss
	 * -128.0950 (28.6% / 4.3%), Met loss -131.0405 (16.7% / 0%), acetyl +42.0106
	 * (28.6% / 14.3%), trimethyl +42.0470 (0% / 0%), carbamylation +43.0058
	 * (15.8% / 10.5%), dioxidation +31.9898 (67.5%, present in only 2 of 9 runs),
	 * trioxidation +47.9847 (11.5%), phospho +79.9663 (0%), propionamide +71.0371
	 * (27.6%), cysteinylation +119.0041 (0%).
	 * <p>
	 * Two further candidates passed every criterion and are still excluded on
	 * principle: -0.9840 (87.9%) sits 19.3 mDa from the zero family's -1 isotope
	 * and the separation is not safe at production tolerances, and an unnamed
	 * +3.9949 (91.1%) would be discovery by the back door. Revisit both only with
	 * a measurement, not an argument.
	 */
	public static final List<ExpectedShift> DEFAULT_SHIFTS = Collections.unmodifiableList(new ArrayList<>(List.of(
			new ExpectedShift("oxidation M/W", 15.99491),
			new ExpectedShift("deamidation N/Q", 0.98402),
			new ExpectedShift("pyroGlu from Q / ammonia loss", -17.02655),
			new ExpectedShift("water loss / pyroGlu from E", -18.01056),
			new ExpectedShift("Na-H", 21.98194),
			new ExpectedShift("Ca-2H", 37.94694),
			new ExpectedShift("carbamidomethyl beyond the fixed mod", 57.02146))));

	/** Thresholds for the four criteria, plus the global cap. */
	public static final class Config {
		/**
		 * Off returns an empty result and changes nothing, so the previous
		 * zero-family-only behaviour stays one flag away.
		 */
		public boolean enabled = true;
		/** Allow-list; defaults to {@link #DEFAULT_SHIFTS}. */
		public List<ExpectedShift> shifts = DEFAULT_SHIFTS;
		/**
		 * Half-width of each shift's window, in Da, around the run's calibrated
		 * position. 8 mDa separates acetyl (+42.0106) from trimethyl (+42.0470),
		 * which are 36 mDa apart and behave completely differently.
		 */
		public double toleranceDa = 0.008;
		/** Criterion 2: minimum target PSMs in the window for it to be testable. */
		public int minTargetsInWindow = 20;
		/** Criterion 3: maximum raw decoys/targets in the window, before any cut. */
		public double maxWindowFdr = 0.50;
		/** Criterion 4: minimum share of accepted peptides also seen unmodified. */
		public double minCorroboration = 0.50;
		/** Global cap on admitted targets, as a fraction of the zero family's. */
		public double maxAdmittedFraction = 0.25;
		/** Per-peak q-value cut; matches {@link UnmodifiedPsmSelector#DEFAULT_FDR}. */
		public double fdrThreshold = UnmodifiedPsmSelector.DEFAULT_FDR;
	}

	/** What happened to one candidate shift, admitted or not. */
	public static final class Verdict {
		public final ExpectedShift shift;
		public final boolean admitted;
		/** Human-readable reason, always populated - "admitted" when it passed. */
		public final String reason;
		public final int targetsInWindow;
		public final int decoysInWindow;
		public final double windowFdr;
		public final int acceptedTargets;
		public final int acceptedDecoys;
		/** Share of accepted target peptides also accepted in the zero family. */
		public final double corroboration;

		Verdict(ExpectedShift shift, boolean admitted, String reason, int targetsInWindow, int decoysInWindow,
				double windowFdr, int acceptedTargets, int acceptedDecoys, double corroboration) {
			this.shift = shift;
			this.admitted = admitted;
			this.reason = reason;
			this.targetsInWindow = targetsInWindow;
			this.decoysInWindow = decoysInWindow;
			this.windowFdr = windowFdr;
			this.acceptedTargets = acceptedTargets;
			this.acceptedDecoys = acceptedDecoys;
			this.corroboration = corroboration;
		}
	}

	public static final class Result {
		/**
		 * PSMs admitted from all accepted shifts, targets and decoys alike. Decoys
		 * are kept for the same reason {@link UnmodifiedPsmSelector} keeps them:
		 * the genome-level FDR downstream needs them.
		 */
		public final List<OpenSearchPsm> admitted;
		/** Every candidate on the allow-list, admitted or rejected, with its reason. */
		public final List<Verdict> verdicts;
		/** Shift name per admitted PSM, by identity; for the audit trail. */
		public final Map<OpenSearchPsm, String> shiftLabels;
		public final int admittedTargets;
		public final int admittedDecoys;

		Result(List<OpenSearchPsm> admitted, List<Verdict> verdicts, Map<OpenSearchPsm, String> shiftLabels,
				int admittedTargets, int admittedDecoys) {
			this.admitted = admitted;
			this.verdicts = verdicts;
			this.shiftLabels = shiftLabels;
			this.admittedTargets = admittedTargets;
			this.admittedDecoys = admittedDecoys;
		}

		/** Shifts that passed all four criteria, in the order they were admitted. */
		public List<ExpectedShift> admittedShifts() {
			List<ExpectedShift> out = new ArrayList<>();
			for (Verdict v : verdicts) {
				if (v.admitted) {
					out.add(v.shift);
				}
			}
			return out;
		}
	}

	private ExpectedPeakSelector() {
	}

	/**
	 * The allow-list entry nearest {@code massDa}, or null when none lies within
	 * {@code toleranceDa} or the mass is not finite.
	 */
	public static ExpectedShift matching(double massDa, List<ExpectedShift> shifts, double toleranceDa) {
		if (!Double.isFinite(massDa) || shifts == null) {
			return null;
		}
		ExpectedShift best = null;
		double bestDistance = Double.MAX_VALUE;
		for (ExpectedShift shift : shifts) {
			double distance = Math.abs(massDa - shift.massDa);
			if (distance <= toleranceDa && distance < bestDistance) {
				bestDistance = distance;
				best = shift;
			}
		}
		return best;
	}

	private static Result empty(List<Verdict> verdicts) {
		return new Result(new ArrayList<>(), verdicts, new IdentityHashMap<>(), 0, 0);
	}

	/**
	 * @param psms     the same rank-1 PSM list {@link UnmodifiedPsmSelector} was
	 *                 given, unmodified. Members of the zero family are skipped
	 *                 here rather than re-tested, so assignment stays exclusive.
	 * @param zero     the zero-family selection for this run; supplies the
	 *                 calibration offset, the membership to exclude, and the
	 *                 accepted peptides that criterion 4 tests against.
	 * @param cfg      thresholds; null uses the defaults.
	 * @param progress optional message sink, mirroring the other stage-1 steps.
	 */
	public static Result select(List<OpenSearchPsm> psms, UnmodifiedPsmSelector.Result zero, Config cfg,
			Consumer<String> progress) {
		if (cfg == null) {
			cfg = new Config();
		}
		List<Verdict> verdicts = new ArrayList<>();
		if (!cfg.enabled || psms == null || psms.isEmpty() || zero == null) {
			return empty(verdicts);
		}
		if (!Double.isFinite(zero.calibrationOffset)) {
			// No zero peak was fitted, so there is no calibrated position to place
			// the windows at. Returning empty leaves the caller with exactly the
			// zero-family behaviour it had before.
			LOGGER.warn("Stage 1 expected peaks: no finite calibration offset from the zero family; skipping.");
			return empty(verdicts);
		}

		// --- exclusive assignment -------------------------------------------
		// The zero family (including the isotope satellites it already folded onto
		// its root) claims its members first. Every remaining PSM goes to its
		// nearest single shift, so no PSM is counted twice and no window can
		// borrow evidence from its neighbour.
		Set<OpenSearchPsm> claimed = Collections.newSetFromMap(new IdentityHashMap<>());
		if (zero.zeroFamily != null) {
			claimed.addAll(zero.zeroFamily);
		}

		List<ExpectedShift> shifts = cfg.shifts == null ? DEFAULT_SHIFTS : cfg.shifts;
		List<List<OpenSearchPsm>> windows = new ArrayList<>();
		for (int i = 0; i < shifts.size(); i++) {
			windows.add(new ArrayList<>());
		}
		for (OpenSearchPsm psm : psms) {
			if (claimed.contains(psm)) {
				continue;
			}
			double observed = psm.getMassDiff() - zero.calibrationOffset;
			int best = -1;
			double bestDistance = Double.MAX_VALUE;
			for (int i = 0; i < shifts.size(); i++) {
				double distance = Math.abs(observed - shifts.get(i).massDa);
				if (distance <= cfg.toleranceDa && distance < bestDistance) {
					bestDistance = distance;
					best = i;
				}
			}
			if (best >= 0) {
				windows.get(best).add(psm);
			}
		}

		// --- test each window, independently --------------------------------
		List<Candidate> passed = new ArrayList<>();
		for (int i = 0; i < shifts.size(); i++) {
			ExpectedShift shift = shifts.get(i);
			List<OpenSearchPsm> window = windows.get(i);
			int targets = 0;
			for (OpenSearchPsm psm : window) {
				if (psm.isTarget()) {
					targets++;
				}
			}
			int decoys = window.size() - targets;

			if (targets < cfg.minTargetsInWindow) {
				verdicts.add(new Verdict(shift, false,
						"peak absent: " + targets + " target PSMs, needs " + cfg.minTargetsInWindow, targets, decoys,
						fdrOf(targets, decoys), 0, 0, Double.NaN));
				continue;
			}
			double windowFdr = fdrOf(targets, decoys);
			if (windowFdr > cfg.maxWindowFdr) {
				verdicts.add(new Verdict(shift, false,
						String.format(Locale.ROOT, "not target-enriched: window decoys/targets %.3f > %.3f", windowFdr,
								cfg.maxWindowFdr),
						targets, decoys, windowFdr, 0, 0, Double.NaN));
				continue;
			}

			// Its own q-value cut, within this peak only.
			double[] scores = new double[window.size()];
			boolean[] isTarget = new boolean[window.size()];
			for (int k = 0; k < window.size(); k++) {
				OpenSearchPsm psm = window.get(k);
				double expect = psm.getExpect();
				scores[k] = Double.isNaN(expect) ? psm.getHyperscore() : -Math.log10(Math.max(expect, 1e-300));
				isTarget[k] = psm.isTarget();
			}
			double[] q = TargetDecoyFdr.qValues(scores, isTarget);

			List<OpenSearchPsm> accepted = new ArrayList<>();
			int acceptedTargets = 0;
			int acceptedDecoys = 0;
			int corroborated = 0;
			int acceptedTargetPeptides = 0;
			// Corroboration is a peptide-level share, so each accepted target
			// peptide votes once however many spectra carry it.
			Set<String> distinct = new HashSet<>();
			for (int k = 0; k < window.size(); k++) {
				if (q[k] > cfg.fdrThreshold) {
					continue;
				}
				OpenSearchPsm psm = window.get(k);
				accepted.add(psm);
				if (isTarget[k]) {
					acceptedTargets++;
					String peptide = psm.getPeptide();
					if (peptide != null && distinct.add(peptide)) {
						acceptedTargetPeptides++;
						if (zero.acceptedTargetPeptides != null && zero.acceptedTargetPeptides.contains(peptide)) {
							corroborated++;
						}
					}
				} else {
					acceptedDecoys++;
				}
			}

			if (acceptedTargets == 0) {
				verdicts.add(new Verdict(shift, false, "no target PSM survived its own q<=" + cfg.fdrThreshold + " cut",
						targets, decoys, windowFdr, 0, acceptedDecoys, Double.NaN));
				continue;
			}
			double corroboration = acceptedTargetPeptides == 0 ? 0.0
					: (double) corroborated / (double) acceptedTargetPeptides;
			if (corroboration < cfg.minCorroboration) {
				verdicts.add(new Verdict(shift, false,
						String.format(Locale.ROOT,
								"uncorroborated: only %.1f%% of its %d accepted peptides were also identified "
										+ "unmodified in this run, needs %.1f%%",
								100.0 * corroboration, acceptedTargetPeptides, 100.0 * cfg.minCorroboration),
						targets, decoys, windowFdr, acceptedTargets, acceptedDecoys, corroboration));
				continue;
			}
			passed.add(new Candidate(shift, accepted, targets, decoys, windowFdr, acceptedTargets, acceptedDecoys,
					corroboration));
		}

		// --- global cap, cleanest peaks first -------------------------------
		// Ordering by window FDR rather than by size means that if the cap does
		// bite, what survives is the best-supported evidence rather than the most
		// abundant.
		passed.sort((a, b) -> Double.compare(a.windowFdr, b.windowFdr));
		int zeroAcceptedTargets = zero.targetsAfterFdr;
		int budget = zeroAcceptedTargets <= 0 ? Integer.MAX_VALUE
				: (int) Math.floor(cfg.maxAdmittedFraction * zeroAcceptedTargets);

		List<OpenSearchPsm> admitted = new ArrayList<>();
		Map<OpenSearchPsm, String> labels = new IdentityHashMap<>();
		int admittedTargets = 0;
		int admittedDecoys = 0;
		for (Candidate candidate : passed) {
			if (admittedTargets + candidate.acceptedTargets > budget) {
				verdicts.add(new Verdict(candidate.shift, false,
						String.format(Locale.ROOT,
								"over the cap: %d admitted targets already, adding %d would exceed %.0f%% of the "
										+ "zero family's %d",
								admittedTargets, candidate.acceptedTargets, 100.0 * cfg.maxAdmittedFraction,
								zeroAcceptedTargets),
						candidate.targets, candidate.decoys, candidate.windowFdr, candidate.acceptedTargets,
						candidate.acceptedDecoys, candidate.corroboration));
				continue;
			}
			for (OpenSearchPsm psm : candidate.accepted) {
				// Normalized onto the zero scale, so a genome-mapping consumer that
				// reads getMassDiff() through OpenSearchMetaPsm sees an
				// unmodified-equivalent peptide, while getMassDiff() on the PSM
				// itself keeps the raw observation for audit.
				psm.setNormalizedMassDiff(psm.getMassDiff() - candidate.shift.massDa);
				psm.setPeakIsotopeOffset(0);
				admitted.add(psm);
				labels.put(psm, candidate.shift.name);
			}
			admittedTargets += candidate.acceptedTargets;
			admittedDecoys += candidate.acceptedDecoys;
			verdicts.add(new Verdict(candidate.shift, true, "admitted", candidate.targets, candidate.decoys,
					candidate.windowFdr, candidate.acceptedTargets, candidate.acceptedDecoys, candidate.corroboration));
		}

		report(verdicts, admittedTargets, admittedDecoys, zeroAcceptedTargets, progress);
		return new Result(admitted, verdicts, labels, admittedTargets, admittedDecoys);
	}

	private static double fdrOf(int targets, int decoys) {
		return targets == 0 ? Double.NaN : (double) decoys / (double) targets;
	}

	private static void report(List<Verdict> verdicts, int admittedTargets, int admittedDecoys,
			int zeroAcceptedTargets, Consumer<String> progress) {
		StringBuilder admittedNames = new StringBuilder();
		int admittedCount = 0;
		for (Verdict v : verdicts) {
			if (!v.admitted) {
				continue;
			}
			if (admittedCount++ > 0) {
				admittedNames.append(", ");
			}
			admittedNames.append(String.format(Locale.ROOT, "%s: %d targets (window FDR %.3f, corroboration %.0f%%)",
					v.shift.toString(), v.acceptedTargets, v.windowFdr, 100.0 * v.corroboration));
		}
		String summary = String.format(Locale.ROOT,
				"Stage 1 expected peaks: %d of %d candidate shifts admitted, adding %d target/%d decoy PSMs "
						+ "(%.1f%% of the zero family's %d accepted targets). %s",
				admittedCount, verdicts.size(), admittedTargets, admittedDecoys,
				zeroAcceptedTargets <= 0 ? 0.0 : 100.0 * admittedTargets / zeroAcceptedTargets, zeroAcceptedTargets,
				admittedCount == 0 ? "None passed." : admittedNames.toString());
		LOGGER.info(summary);
		if (progress != null) {
			progress.accept(summary);
		}
		for (Verdict v : verdicts) {
			if (v.admitted) {
				continue;
			}
			LOGGER.debug("Stage 1 expected peaks: rejected " + v.shift + " - " + v.reason);
		}
	}

	/** A shift that passed all four criteria, held until the cap is applied. */
	private static final class Candidate {
		final ExpectedShift shift;
		final List<OpenSearchPsm> accepted;
		final int targets;
		final int decoys;
		final double windowFdr;
		final int acceptedTargets;
		final int acceptedDecoys;
		final double corroboration;

		Candidate(ExpectedShift shift, List<OpenSearchPsm> accepted, int targets, int decoys, double windowFdr,
				int acceptedTargets, int acceptedDecoys, double corroboration) {
			this.shift = shift;
			this.accepted = accepted;
			this.targets = targets;
			this.decoys = decoys;
			this.windowFdr = windowFdr;
			this.acceptedTargets = acceptedTargets;
			this.acceptedDecoys = acceptedDecoys;
			this.corroboration = corroboration;
		}
	}
}
