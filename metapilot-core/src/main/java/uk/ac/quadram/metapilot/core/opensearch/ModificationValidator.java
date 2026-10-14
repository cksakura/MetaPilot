package uk.ac.quadram.metapilot.core.opensearch;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.Set;
import java.util.function.Consumer;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import uk.ac.quadram.metapilot.core.opensearch.ExpectedPeakSelector.ExpectedShift;

/**
 * Per-modification validation of scored open-search PSMs.
 * <p>
 * Every PSM reaching this class is routed into exactly one error-controlled
 * population. Nothing is discarded without a target-decoy test:
 * <ol>
 * <li><b>peak-validated</b> — the PSM lies in a fitted delta-mass peak that the
 * screen accepted. Cut by a local FDR walk <em>within that peak</em>
 * ({@link Config#localFdr}), and reported as evidence for the modification the
 * peak represents.</li>
 * <li><b>rejected-peak identifications</b> — the PSM lies in a fitted peak the
 * screen did not accept. Cut by its own identification FDR walk
 * ({@link Config#unexplainedFdr}).</li>
 * <li><b>no-peak identifications</b> — the delta mass falls outside every
 * fitted peak ({@code peakId == -1}). Cut by a separate walk at the same
 * threshold.</li>
 * </ol>
 *
 * <h2>Why a rejected peak demotes rather than deletes</h2>
 * Rejecting a peak is a statement about the <em>mass shift</em> — "this offset
 * is not a modification I can support" — not about the spectra that carry it.
 * Deleting rejected-peak PSMs would make peak membership decide
 * identification, and would do so perversely: of two PSMs with the same
 * peptide and the same score, the one whose mass clustered with others would
 * be destroyed while the one whose mass fell in the background survived.
 *
 * <h2>Why the two identification-only populations are walked separately</h2>
 * They are exported under different scope labels ({@code REJECTED_PEAK_GROUP},
 * {@code NO_PEAK_GROUP}), and a q-value is only a statement about the
 * population it was computed in. Rejected-peak PSMs include the bulk of the
 * clean, high-scoring members of peaks that make no modification claim, so a
 * pooled walk lets them dilute the decoys of the background population and
 * reports the background at an FDR it does not have.
 *
 * <h2>Why the peak screen is not redundant with the walk</h2>
 * Without it a pure-noise peak still contributes its luckiest PSMs: if its top
 * two happen to be targets then {@code decoys/targets = 0/2 = 0} and they pass
 * as a validated modification. The screen's job is to stop that claim being
 * made, not to decide which spectra are identified.
 *
 * <h2>Calibration caveat</h2>
 * Peaks are <em>defined</em> by fitting a target-dominated histogram, so
 * conditioning on "is a fitted peak" already pre-selects for target enrichment
 * and the screen's decoy fractions are optimistic rather than strictly
 * calibrated. That is acceptable because the screen is not the error control:
 * FDR control comes from the walks, which use the 1:1 decoys directly.
 */
public final class ModificationValidator {

	private static final Logger LOGGER = LogManager.getLogger(ModificationValidator.class);

	/** Two-sided 95% normal quantile for the Wilson bound. */
	private static final double WILSON_Z = 1.959964;

	/**
	 * Background decoy fraction used when no other non-zero peak exists to
	 * measure one: under a 1:1 target-decoy search random matches split evenly.
	 */
	static final double NULL_DECOY_FRACTION = 0.5;

	public static final class Config {
		/**
		 * Half-width, in Da, of the window around zero in which the peak nearest
		 * zero is taken as the unmodified family. Matches
		 * {@link DeltaMassEvidenceGate.Config#zeroWindowDa}.
		 */
		public double zeroWindowDa = 0.05;
		/** Masses that may be accepted as modifications; everything else is unexplained. */
		public List<ExpectedShift> knownShifts = ExpectedPeakSelector.DEFAULT_SHIFTS;
		/** Matching tolerance, in Da, between a peak centroid and a known shift. */
		public double knownShiftToleranceDa = 0.008;
		/** Local FDR threshold within each accepted peak. */
		public double localFdr = 0.01;
		/** Identification FDR threshold for each identification-only population. */
		public double unexplainedFdr = 0.01;
	}

	/** Which hypothesis a peak's mass represents. */
	public enum MassClass {
		/** The peak nearest zero: the unmodified hypothesis. */
		ZERO,
		/** Within tolerance of an allow-list shift. */
		KNOWN,
		/** Any other mass. */
		UNEXPLAINED
	}

	public static final class PeakVerdict {
		public final int peakId;
		public final double centroid;
		public final MassClass massClass;
		/** Allow-list entry name for a {@link MassClass#KNOWN} peak, otherwise empty. */
		public final String knownShift;
		public final int targets;
		public final int decoys;
		/** The peak's own decoy fraction among rank-1 in-peak PSMs. */
		public final double decoyFraction;
		/** Wilson 95% upper bound on {@link #decoyFraction}. */
		public final double decoyUpperBound;
		/**
		 * Decoy fraction of the rank-1 PSMs in every other non-zero peak - the
		 * reference a known shift must beat. NaN for the zero family.
		 */
		public final double backgroundDecoyFraction;
		/**
		 * True when {@link DeltaMassEvidenceGate} scored enough members and most
		 * gained from the shifted fragment ions. Recorded, not used to accept.
		 */
		public final boolean fragmentSupported;
		public final boolean accepted;
		public final String reason;
		public final int psmsRetained;

		PeakVerdict(int peakId, double centroid, MassClass massClass, String knownShift, int targets, int decoys,
				double decoyFraction, double decoyUpperBound, double backgroundDecoyFraction,
				boolean fragmentSupported, boolean accepted, String reason, int psmsRetained) {
			this.peakId = peakId;
			this.centroid = centroid;
			this.massClass = massClass;
			this.knownShift = knownShift;
			this.targets = targets;
			this.decoys = decoys;
			this.decoyFraction = decoyFraction;
			this.decoyUpperBound = decoyUpperBound;
			this.backgroundDecoyFraction = backgroundDecoyFraction;
			this.fragmentSupported = fragmentSupported;
			this.accepted = accepted;
			this.reason = reason;
			this.psmsRetained = psmsRetained;
		}

		PeakVerdict withdrawn(String why) {
			return new PeakVerdict(peakId, centroid, massClass, knownShift, targets, decoys, decoyFraction,
					decoyUpperBound, backgroundDecoyFraction, fragmentSupported, false, why, 0);
		}

		PeakVerdict withRetained(int retained) {
			return new PeakVerdict(peakId, centroid, massClass, knownShift, targets, decoys, decoyFraction,
					decoyUpperBound, backgroundDecoyFraction, fragmentSupported, accepted, reason, retained);
		}
	}

	public static final class Result {
		/** Union used for peptide/protein/genome reporting. */
		public final List<OpenSearchPsm> validated;
		/** PSMs supported by an accepted fitted delta-mass hypothesis. */
		public final List<OpenSearchPsm> peakValidated;
		/**
		 * Confident peptide IDs with no accepted mass explanation - outside every
		 * fitted peak, or inside a peak the screen rejected, each population cut
		 * by its own walk. Identifications, not modification claims.
		 */
		public final List<OpenSearchPsm> unexplainedIdentifications;
		public final List<PeakVerdict> verdicts;
		/**
		 * Open PSMs the variant outputs read: accepted plus unexplained. The same list as
		 * {@link #validated} unless the closed stage took the unexplained tier out of it.
		 */
		public final List<OpenSearchPsm> diagnostic;

		Result(List<OpenSearchPsm> validated, List<OpenSearchPsm> peakValidated,
				List<OpenSearchPsm> unexplainedIdentifications, List<PeakVerdict> verdicts) {
			this(validated, peakValidated, unexplainedIdentifications, verdicts, validated);
		}

		public Result(List<OpenSearchPsm> validated, List<OpenSearchPsm> peakValidated,
				List<OpenSearchPsm> unexplainedIdentifications, List<PeakVerdict> verdicts,
				List<OpenSearchPsm> diagnostic) {
			this.validated = validated;
			this.peakValidated = peakValidated;
			this.unexplainedIdentifications = unexplainedIdentifications;
			this.verdicts = verdicts;
			this.diagnostic = diagnostic;
		}

		public int acceptedPeakCount() {
			int n = 0;
			for (PeakVerdict v : verdicts) {
				if (v.accepted) {
					n++;
				}
			}
			return n;
		}
	}

	private ModificationValidator() {
	}

	/** What {@link #filterPeaks} decided, before any PSM is scored. */
	public static final class Gate {
		public final Set<Integer> acceptedPeakIds;
		public final List<PeakVerdict> verdicts;
		/** Files whose PSMs were used to test; empty when the gate is not independent. */
		public final Set<String> gateFiles;

		Gate(Set<Integer> acceptedPeakIds, List<PeakVerdict> verdicts, Set<String> gateFiles) {
			this.acceptedPeakIds = acceptedPeakIds;
			this.verdicts = verdicts;
			this.gateFiles = gateFiles;
		}

		public int acceptedCount() {
			int n = 0;
			for (PeakVerdict v : verdicts) {
				if (v.accepted) {
					n++;
				}
			}
			return n;
		}

		/**
		 * The same gate with some peaks withdrawn from the accepted set.
		 * <p>
		 * Used for peaks that pass the screen but are not modification hypotheses
		 * at all - a peak whose delta is explained by a single amino-acid
		 * substitution rather than a mass change on the assigned peptide. Their
		 * PSMs fall back to identification-only, which is exactly what withdrawing
		 * them from {@code acceptedPeakIds} does. The withdrawal is recorded in
		 * the verdict's reason.
		 */
		public Gate withoutPeaks(Set<Integer> withdrawn) {
			if (withdrawn == null || withdrawn.isEmpty()) {
				return this;
			}
			Set<Integer> kept = new HashSet<>(acceptedPeakIds);
			kept.removeAll(withdrawn);
			List<PeakVerdict> revised = new ArrayList<>(verdicts.size());
			for (PeakVerdict v : verdicts) {
				revised.add(withdrawn.contains(v.peakId) && v.accepted
						? v.withdrawn("withdrawn: explained by an amino-acid substitution")
						: v);
			}
			return new Gate(kept, revised, gateFiles);
		}
	}

	/**
	 * Decides which fitted peaks may be reported as modifications.
	 *
	 * <h2>The rule</h2>
	 * <ul>
	 * <li>The peak nearest zero (within {@link Config#zeroWindowDa}) is the
	 * unmodified hypothesis and is always accepted.</li>
	 * <li>A peak within {@link Config#knownShiftToleranceDa} of an allow-list
	 * shift is accepted when the Wilson 95% upper bound on its rank-1 decoy
	 * fraction lies below the decoy fraction of all other non-zero peaks
	 * combined: the shift must be target-enriched against the in-peak
	 * background, with the bound making a small peak earn it.</li>
	 * <li>Any other mass is not accepted. Its PSMs stay identifications and no
	 * modification claim is made.</li>
	 * </ul>
	 *
	 * <h2>Why not a cumulative walk to a target FDR</h2>
	 * The walk this replaced admitted peaks in order of decoy fraction until the
	 * admitted pool reached 1.5%. On PXD050781 it rejected the zero peak itself -
	 * 1.93% decoys over 592,038 rank-1 PSMs for MSFragger, 1.59% for MSFragger +
	 * Sage - because hundreds of small unexplained peaks with 0 decoys out of a
	 * few dozen PSMs were cleaner and filled the budget first. Oxidation was
	 * never reached. A pooled walk lets large clean peaks carry small dirty ones
	 * and stops on a single large one, so its verdicts depend on what else was
	 * fitted rather than on the peak.
	 *
	 * <h2>Why the background, not a multiple of the zero peak</h2>
	 * Decoy fraction measures identification quality, not whether a mass is a
	 * modification. On PXD050781 the upper bounds of real known shifts ranged
	 * from 0.7x to 5.1x the zero peak's decoy fraction with no gap separating
	 * them from junk, so any fixed multiple is arbitrary. The non-zero in-peak
	 * background (about 19% decoys on both workflows) is measured in each run
	 * and carries over between datasets; it is computed leave-one-out so a large
	 * peak is not held to a reference it dominates.
	 *
	 * <h2>Why unexplained masses make no claim</h2>
	 * Target enrichment cannot distinguish an unknown modification from an
	 * amino-acid substitution, a chimeric spectrum or a wrong peptide carrying a
	 * consistent offset, and neither can fragment gain: on PXD050781 known shifts
	 * and substitution-shaped peaks gain from shifted fragments at similar rates.
	 * Until a test exists that separates them, the PSMs are reported as
	 * identifications and the fragment-support status is recorded for audit.
	 *
	 * <h2>Why rank 1 only</h2>
	 * Nothing has been rescored yet, so the engine's own ranking is the only
	 * ordering available. Lower ranks are decoy-enriched everywhere, so counting
	 * them would measure rank depth rather than peak quality.
	 *
	 * @param fragmentSupported peak ids {@link DeltaMassEvidenceGate} found
	 *                          supported by shifted fragment ions; recorded in
	 *                          the verdicts. Null means not measured.
	 */
	public static Gate filterPeaks(List<OpenSearchPsm> psms, Config cfg, Set<Integer> fragmentSupported,
			Consumer<String> progress) {
		if (cfg == null) {
			cfg = new Config();
		}
		Map<Integer, int[]> counts = new LinkedHashMap<>();
		Map<Integer, Double> centroids = new HashMap<>();
		int used = 0;
		for (OpenSearchPsm psm : psms) {
			if (psm.getPeakId() < 0 || psm.getHitRank() != 1) {
				continue;
			}
			used++;
			int[] tallies = counts.computeIfAbsent(psm.getPeakId(), k -> new int[2]);
			tallies[psm.isTarget() ? 0 : 1]++;
			centroids.putIfAbsent(psm.getPeakId(), psm.getPeakCentroid());
		}

		List<PeakVerdict> verdicts = screen(counts, centroids, cfg, fragmentSupported);
		Set<Integer> accepted = new LinkedHashSet<>();
		int known = 0;
		int unexplained = 0;
		int unexplainedSupported = 0;
		for (PeakVerdict v : verdicts) {
			if (v.accepted) {
				accepted.add(v.peakId);
			}
			if (v.massClass == MassClass.KNOWN) {
				known++;
			} else if (v.massClass == MassClass.UNEXPLAINED) {
				unexplained++;
				if (v.fragmentSupported) {
					unexplainedSupported++;
				}
			}
		}
		boolean zeroFound = false;
		for (PeakVerdict v : verdicts) {
			zeroFound |= v.massClass == MassClass.ZERO;
		}
		String note = String.format(Locale.ROOT,
				"Peak screen: %d/%d peaks accepted from %,d rank-1 in-peak PSMs (zero family %s; %d/%d known "
						+ "shifts target-enriched; %d unexplained masses identification-only, %d of them "
						+ "fragment-supported)",
				accepted.size(), verdicts.size(), used, zeroFound ? "accepted" : "not found",
				accepted.size() - (zeroFound ? 1 : 0), known, unexplained, unexplainedSupported);
		LOGGER.info(note);
		if (progress != null) {
			progress.accept(note);
		}
		return new Gate(accepted, verdicts, Collections.emptySet());
	}

	/**
	 * Applies the rule documented on {@link #filterPeaks} to per-peak target and
	 * decoy counts. Shared with {@link #validate}'s no-gate path so the two
	 * cannot diverge.
	 */
	static List<PeakVerdict> screen(Map<Integer, int[]> counts, Map<Integer, Double> centroids, Config cfg,
			Set<Integer> fragmentSupported) {
		int zeroPeakId = zeroPeakId(centroids, cfg.zeroWindowDa);
		long nonZeroTargets = 0;
		long nonZeroDecoys = 0;
		for (Map.Entry<Integer, int[]> e : counts.entrySet()) {
			if (e.getKey() != zeroPeakId) {
				nonZeroTargets += e.getValue()[0];
				nonZeroDecoys += e.getValue()[1];
			}
		}

		List<PeakVerdict> verdicts = new ArrayList<>(counts.size());
		for (Map.Entry<Integer, int[]> e : counts.entrySet()) {
			int peakId = e.getKey();
			int targets = e.getValue()[0];
			int decoys = e.getValue()[1];
			int n = targets + decoys;
			double centroid = centroids.getOrDefault(peakId, Double.NaN);
			double fraction = n == 0 ? 1.0 : decoys / (double) n;
			double upper = wilsonUpperBound(decoys, n);
			boolean supported = fragmentSupported != null && fragmentSupported.contains(peakId);

			if (peakId == zeroPeakId) {
				verdicts.add(new PeakVerdict(peakId, centroid, MassClass.ZERO, "", targets, decoys, fraction, upper,
						Double.NaN, supported, true, "zero family: the unmodified hypothesis", 0));
				continue;
			}

			long otherTargets = nonZeroTargets - targets;
			long otherDecoys = nonZeroDecoys - decoys;
			double background = otherTargets + otherDecoys == 0 ? NULL_DECOY_FRACTION
					: otherDecoys / (double) (otherTargets + otherDecoys);
			ExpectedShift shift = ExpectedPeakSelector.matching(centroid, cfg.knownShifts, cfg.knownShiftToleranceDa);
			if (shift == null) {
				verdicts.add(new PeakVerdict(peakId, centroid, MassClass.UNEXPLAINED, "", targets, decoys, fraction,
						upper, background, supported, false, "unexplained mass: identification only", 0));
				continue;
			}
			boolean enriched = upper < background;
			verdicts.add(new PeakVerdict(peakId, centroid, MassClass.KNOWN, shift.name, targets, decoys, fraction,
					upper, background, supported, enriched,
					enriched ? "known shift, target-enriched against the background"
							: "known shift, not target-enriched against the background",
					0));
		}
		verdicts.sort((a, b) -> Integer.compare(a.peakId, b.peakId));
		return verdicts;
	}

	/** The peak with the smallest finite |centroid| inside {@code windowDa}, or -1. */
	static int zeroPeakId(Map<Integer, Double> centroids, double windowDa) {
		int best = -1;
		double bestDistance = Double.MAX_VALUE;
		for (Map.Entry<Integer, Double> e : centroids.entrySet()) {
			double c = e.getValue() == null ? Double.NaN : e.getValue();
			if (!Double.isFinite(c)) {
				continue;
			}
			double distance = Math.abs(c);
			if (distance < windowDa && distance < bestDistance) {
				bestDistance = distance;
				best = e.getKey();
			}
		}
		return best;
	}

	/** Wilson score 95% upper bound on a binomial proportion; 1.0 when {@code n == 0}. */
	static double wilsonUpperBound(int successes, int n) {
		if (n <= 0) {
			return 1.0;
		}
		double p = successes / (double) n;
		double z2 = WILSON_Z * WILSON_Z;
		double centre = p + z2 / (2.0 * n);
		double half = WILSON_Z * Math.sqrt(p * (1 - p) / n + z2 / (4.0 * n * n));
		return Math.min(1.0, (centre + half) / (1 + z2 / n));
	}

	public static Result validate(List<OpenSearchPsm> psms, Config cfg) {
		return validate(psms, cfg, null);
	}

	/**
	 * <b>Contract:</b> {@code psms} must already contain at most one PSM per
	 * spectrum — see {@link SpectrumLevelCollapse}. Both the peak screen and the
	 * per-group FDR walks are target-decoy statistics, and a target-decoy count is
	 * only interpretable when each row is one independent identification. Passing
	 * the raw multi-candidate list is silently wrong rather than an error: one
	 * spectrum can be accepted twice, and because a correct spectrum's runner-up
	 * candidates are correlated targets while decoys stay independent, the true
	 * FDR runs above the reported one.
	 *
	 * @param psms PSMs carrying a classifier score ({@code classLogit}) and a
	 *             delta-mass peak assignment ({@code peakId})
	 * @param gate verdicts from {@link #filterPeaks}, decided before scoring. When
	 *             supplied, peaks are not re-screened on the classifier-collapsed
	 *             set. Null screens the PSMs in hand, so a missing gate never
	 *             silently disables filtering.
	 */
	public static Result validate(List<OpenSearchPsm> psms, Config cfg, Gate gate) {
		if (cfg == null) {
			cfg = new Config();
		}

		Map<Integer, List<OpenSearchPsm>> byPeak = new LinkedHashMap<>();
		List<OpenSearchPsm> noPeak = new ArrayList<>();
		for (OpenSearchPsm psm : psms) {
			if (psm.getPeakId() < 0) {
				noPeak.add(psm);
			} else {
				byPeak.computeIfAbsent(psm.getPeakId(), k -> new ArrayList<>()).add(psm);
			}
		}

		if (gate == null) {
			Map<Integer, int[]> localCounts = new LinkedHashMap<>();
			Map<Integer, Double> localCentroids = new HashMap<>();
			for (Map.Entry<Integer, List<OpenSearchPsm>> e : byPeak.entrySet()) {
				int[] tallies = new int[2];
				for (OpenSearchPsm psm : e.getValue()) {
					tallies[psm.isTarget() ? 0 : 1]++;
					localCentroids.putIfAbsent(e.getKey(), psm.getPeakCentroid());
				}
				localCounts.put(e.getKey(), tallies);
			}
			List<PeakVerdict> local = screen(localCounts, localCentroids, cfg, null);
			Set<Integer> accepted = new LinkedHashSet<>();
			for (PeakVerdict v : local) {
				if (v.accepted) {
					accepted.add(v.peakId);
				}
			}
			gate = new Gate(accepted, local, Collections.emptySet());
		}
		Map<Integer, PeakVerdict> gateVerdicts = new HashMap<>();
		for (PeakVerdict v : gate.verdicts) {
			gateVerdicts.put(v.peakId, v);
		}

		List<Integer> peakIds = new ArrayList<>(byPeak.keySet());
		peakIds.sort(null);

		List<OpenSearchPsm> peakValidated = new ArrayList<>();
		List<OpenSearchPsm> rejectedPeakPsms = new ArrayList<>();
		List<PeakVerdict> verdicts = new ArrayList<>();
		for (int peakId : peakIds) {
			List<OpenSearchPsm> group = byPeak.get(peakId);
			boolean accepted = gate.acceptedPeakIds.contains(peakId);
			// Report the numbers the decision was made on - the gate's rank-1
			// counts, not this collapsed set's - or the audit file would contradict
			// the verdict beside it.
			PeakVerdict verdict = gateVerdicts.get(peakId);
			if (verdict == null) {
				verdict = untested(peakId, group, cfg);
			}
			int retained = 0;
			if (accepted) {
				List<OpenSearchPsm> kept = applyLocalFdr(group, cfg.localFdr);
				peakValidated.addAll(kept);
				retained = kept.size();
			} else {
				rejectedPeakPsms.addAll(group);
			}
			verdicts.add(verdict.withRetained(retained));
		}

		// Identification-level walks, one per exported population. Neither is a
		// peak gate nor evidence that the masses within it share one cause.
		List<OpenSearchPsm> rejectedPeakIds = applyLocalFdr(rejectedPeakPsms, cfg.unexplainedFdr);
		List<OpenSearchPsm> noPeakIds = applyLocalFdr(noPeak, cfg.unexplainedFdr);
		List<OpenSearchPsm> unexplainedIdentifications = new ArrayList<>(rejectedPeakIds.size() + noPeakIds.size());
		unexplainedIdentifications.addAll(rejectedPeakIds);
		unexplainedIdentifications.addAll(noPeakIds);

		List<OpenSearchPsm> validated = new ArrayList<>(peakValidated.size() + unexplainedIdentifications.size());
		validated.addAll(peakValidated);
		validated.addAll(unexplainedIdentifications);

		int acceptedPeaks = 0;
		for (PeakVerdict v : verdicts) {
			if (v.accepted) {
				acceptedPeaks++;
			}
		}
		LOGGER.info(String.format(Locale.ROOT,
				"Modification validation: %d/%d fitted peaks accepted; %,d PSMs retained at %.2f%% local FDR; "
						+ "identification-only at %.2f%% FDR per population: %,d/%,d from %d rejected peaks, "
						+ "%,d/%,d outside every peak",
				acceptedPeaks, verdicts.size(), peakValidated.size(), cfg.localFdr * 100, cfg.unexplainedFdr * 100,
				rejectedPeakIds.size(), rejectedPeakPsms.size(), verdicts.size() - acceptedPeaks, noPeakIds.size(),
				noPeak.size()));

		return new Result(validated, peakValidated, unexplainedIdentifications, verdicts);
	}

	/**
	 * Open stage of the closed-then-open identification
	 * ({@code docs/developer/closed-then-open-identification-plan.md}): {@code remaining} holds the
	 * collapsed open winners of the spectra the closed stage did not accept. Two pooled walks, one
	 * over the PSMs of the accepted zero peak (unmodified) and one over those of accepted
	 * allow-list peaks (known shift); their accepted PSMs are the identifications. Everything else
	 * is an unexplained shift: walked as {@link #validate} does, for {@code unexplained_psms.tsv}
	 * and the variant outputs, but not identified.
	 * <p>
	 * Pooled rather than per peak because the closed stage has already taken most confident
	 * spectra: a per-peak walk over what is left is a few decoys on small counts.
	 */
	public static Result validateRemaining(List<OpenSearchPsm> remaining, Config cfg, Gate gate) {
		if (cfg == null) {
			cfg = new Config();
		}
		Map<Integer, PeakVerdict> gateVerdicts = new HashMap<>();
		for (PeakVerdict v : gate.verdicts) {
			gateVerdicts.put(v.peakId, v);
		}
		Map<Integer, List<OpenSearchPsm>> byPeak = new LinkedHashMap<>();
		List<OpenSearchPsm> unmodified = new ArrayList<>();
		List<OpenSearchPsm> knownShift = new ArrayList<>();
		List<OpenSearchPsm> rejectedPeakPsms = new ArrayList<>();
		List<OpenSearchPsm> noPeak = new ArrayList<>();
		for (OpenSearchPsm psm : remaining) {
			int peakId = psm.getPeakId();
			if (peakId < 0) {
				noPeak.add(psm);
				continue;
			}
			byPeak.computeIfAbsent(peakId, k -> new ArrayList<>()).add(psm);
			PeakVerdict v = gateVerdicts.get(peakId);
			if (!gate.acceptedPeakIds.contains(peakId) || v == null) {
				rejectedPeakPsms.add(psm);
			} else if (v.massClass == MassClass.ZERO) {
				unmodified.add(psm);
			} else if (v.massClass == MassClass.KNOWN) {
				knownShift.add(psm);
			} else {
				rejectedPeakPsms.add(psm);
			}
		}
		List<OpenSearchPsm> keptUnmodified = applyLocalFdr(unmodified, cfg.localFdr);
		List<OpenSearchPsm> keptKnown = applyLocalFdr(knownShift, cfg.localFdr);
		List<OpenSearchPsm> peakValidated = new ArrayList<>(keptUnmodified.size() + keptKnown.size());
		peakValidated.addAll(keptUnmodified);
		peakValidated.addAll(keptKnown);

		Map<Integer, Integer> retained = new HashMap<>();
		for (OpenSearchPsm psm : peakValidated) {
			retained.merge(psm.getPeakId(), 1, Integer::sum);
		}
		List<Integer> peakIds = new ArrayList<>(byPeak.keySet());
		peakIds.sort(null);
		List<PeakVerdict> verdicts = new ArrayList<>();
		for (int peakId : peakIds) {
			PeakVerdict verdict = gateVerdicts.get(peakId);
			if (verdict == null) {
				verdict = untested(peakId, byPeak.get(peakId), cfg);
			}
			verdicts.add(verdict.withRetained(retained.getOrDefault(peakId, 0)));
		}

		List<OpenSearchPsm> rejectedPeakIds = applyLocalFdr(rejectedPeakPsms, cfg.unexplainedFdr);
		List<OpenSearchPsm> noPeakIds = applyLocalFdr(noPeak, cfg.unexplainedFdr);
		List<OpenSearchPsm> unexplained = new ArrayList<>(rejectedPeakIds.size() + noPeakIds.size());
		unexplained.addAll(rejectedPeakIds);
		unexplained.addAll(noPeakIds);
		List<OpenSearchPsm> diagnostic = new ArrayList<>(peakValidated.size() + unexplained.size());
		diagnostic.addAll(peakValidated);
		diagnostic.addAll(unexplained);

		LOGGER.info(String.format(Locale.ROOT,
				"Open stage over %,d spectra the closed stage left: %,d/%,d unmodified and %,d/%,d known-shift "
						+ "PSMs accepted at %.2f%% FDR; %,d unexplained-shift PSMs pass their walks and are not "
						+ "identifications",
				remaining.size(), keptUnmodified.size(), unmodified.size(), keptKnown.size(), knownShift.size(),
				cfg.localFdr * 100, unexplained.size()));
		return new Result(new ArrayList<>(peakValidated), peakValidated, unexplained, verdicts, diagnostic);
	}

	/** Mass class of each peak the gate decided, by peak id. */
	public static Map<Integer, MassClass> massClasses(Gate gate) {
		Map<Integer, MassClass> classes = new HashMap<>();
		for (PeakVerdict v : gate.verdicts) {
			classes.put(v.peakId, v.massClass);
		}
		return classes;
	}

	/**
	 * Verdict for a peak present after the spectrum collapse but absent from the
	 * gate - none of its rank-1 PSMs reached the screen, so nothing was decided
	 * and it is not accepted.
	 */
	private static PeakVerdict untested(int peakId, List<OpenSearchPsm> group, Config cfg) {
		int targets = 0;
		int decoys = 0;
		double centroid = Double.NaN;
		for (OpenSearchPsm psm : group) {
			if (psm.isTarget()) {
				targets++;
			} else {
				decoys++;
			}
			if (!Double.isFinite(centroid)) {
				centroid = psm.getPeakCentroid();
			}
		}
		int n = targets + decoys;
		ExpectedShift shift = ExpectedPeakSelector.matching(centroid, cfg.knownShifts, cfg.knownShiftToleranceDa);
		return new PeakVerdict(peakId, centroid, shift == null ? MassClass.UNEXPLAINED : MassClass.KNOWN,
				shift == null ? "" : shift.name, targets, decoys, n == 0 ? 1.0 : decoys / (double) n,
				wilsonUpperBound(decoys, n), Double.NaN, false, false, "not screened: no rank-1 PSMs at the gate", 0);
	}

	/**
	 * Cuts a PSM list at the deepest rank whose monotonised q-value still
	 * satisfies the threshold — not at the first rank to breach it, which a
	 * single unlucky decoy near the top can truncate prematurely.
	 */
	private static List<OpenSearchPsm> applyLocalFdr(List<OpenSearchPsm> group, double localFdr) {
		int n = group.size();
		double[] scores = new double[n];
		boolean[] isTarget = new boolean[n];
		for (int i = 0; i < n; i++) {
			OpenSearchPsm psm = group.get(i);
			scores[i] = psm.getComparableScore();
			isTarget[i] = psm.isTarget();
		}

		double[] q = TargetDecoyFdr.qValues(scores, isTarget);
		List<OpenSearchPsm> kept = new ArrayList<>();
		for (int i = 0; i < n; i++) {
			group.get(i).setQValue(q[i]);
			if (q[i] <= localFdr) {
				kept.add(group.get(i));
			}
		}
		return kept;
	}
}
