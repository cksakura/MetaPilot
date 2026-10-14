package uk.ac.quadram.metapilot.core.opensearch;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

/**
 * Non-parametric delta-mass mode discovery used by {@link DeltaMassPeakFitter}.
 *
 * <p>The implementation operates on the binned distribution because a direct
 * KDE over hundreds of thousands of PSMs and a 650 Da axis is unnecessarily
 * expensive. Convolution with a Gaussian kernel is mathematically equivalent
 * to evaluating a KDE on the histogram grid. Several bandwidths are evaluated;
 * a candidate must persist across scales and, when file labels are supplied,
 * recur in independent files.</p>
 *
 * <p>KDE discovers modes and empirical valleys only. It does not validate a
 * modification: family-level target/decoy and localization validation remain
 * downstream responsibilities.</p>
 */
final class DeltaMassKdePeakFinder {
	private static final Logger LOGGER = LogManager.getLogger(DeltaMassKdePeakFinder.class);

	private DeltaMassKdePeakFinder() {
	}

	static List<DeltaMassPeakFitter.Peak> discover(double[] masses, boolean[] discoveryMask,
			double[] discoveryWeights, boolean[] isTarget, String[] fileGroups,
			DeltaMassPeakFitter.Config cfg) {
		return discover(masses, discoveryMask, discoveryWeights, isTarget, fileGroups, cfg, null);
	}

	static List<DeltaMassPeakFitter.Peak> discover(double[] masses, boolean[] discoveryMask,
			double[] discoveryWeights, boolean[] isTarget, String[] fileGroups,
			DeltaMassPeakFitter.Config cfg, java.util.function.Consumer<String> progress) {
		int bins = (int) Math.ceil((cfg.windowUpper - cfg.windowLower) / cfg.binWidth);
		double[] target = new double[bins];
		double[] decoy = new double[bins];
		Map<String, double[]> targetByFile = new LinkedHashMap<>();
		Map<String, double[]> decoyByFile = new LinkedHashMap<>();

		for (int i = 0; i < masses.length; i++) {
			int bin = binOf(masses[i], cfg);
			if (bin < 0 || bin >= bins) {
				continue;
			}
			boolean targetRow = isTarget == null || isTarget[i];
			boolean discover = discoveryMask == null || discoveryMask[i];
			if (targetRow && !discover) {
				continue;
			}
			// detectionRows gates BOTH classes, unlike discoveryMask above which
			// gates targets only. It carries rank-1 membership: before rescoring the
			// engine's own ranking is the only ordering available, so admitting
			// ranks 2-5 here would contradict the evidence the step relies on.
			// Measured on 12mix, 57.4% of rank-1 candidates fall inside a fitted
			// peak against 16.9% of ranks 2-5, and under per-spectrum weighting
			// those lower ranks carried ~80% of each spectrum's vote into diffuse
			// background. Excluding only their targets would bias the
			// target-minus-decoy excess the modes are found in, so both go.
			if (cfg.detectionRows != null && !cfg.detectionRows[i]) {
				continue;
			}
			double weight = discoveryWeights == null ? 1.0 : discoveryWeights[i];
			if (!Double.isFinite(weight) || weight < 0.0) {
				throw new IllegalArgumentException("discoveryWeights must be finite and non-negative");
			}
			// Existing callers assign zero weight to decoys because the old fitter used
			// a targets-only histogram. Decoys need one spectrum-equivalent vote for
			// target-excess KDE unless an explicit positive weight was provided.
			if (!targetRow && weight == 0.0) {
				weight = 1.0;
			}
			(targetRow ? target : decoy)[bin] += weight;

			if (fileGroups != null && fileGroups[i] != null) {
				Map<String, double[]> map = targetRow ? targetByFile : decoyByFile;
				map.computeIfAbsent(fileGroups[i], key -> new double[bins])[bin] += weight;
			}
		}

		double[] bandwidths = cfg.kdeBandwidthsDa();
		List<double[]> densities = new ArrayList<>(bandwidths.length);
		List<int[]> modesByScale = new ArrayList<>(bandwidths.length);
		double[] decoyBase = null;
		for (int scale = 0; scale < bandwidths.length; scale++) {
			double[] targetDensity = gaussianSmooth(target, bandwidths[scale] / cfg.binWidth);
			double[] decoyDensity = gaussianSmooth(decoy, bandwidths[scale] / cfg.binWidth);
			double[] excess = new double[bins];
			for (int bin = 0; bin < bins; bin++) {
				excess[bin] = Math.max(0.0, targetDensity[bin] - decoyDensity[bin]);
			}
			densities.add(excess);
			modesByScale.add(localMaxima(excess));
			if (scale == 0) {
				decoyBase = decoyDensity;
			}
		}

		double[] base = densities.get(0);
		List<Candidate> candidates = new ArrayList<>();
		for (int apex : modesByScale.get(0)) {
			int persistent = 1;
			for (int scale = 1; scale < bandwidths.length; scale++) {
				double toleranceDa = Math.max(2.0 * cfg.binWidth, 1.5 * bandwidths[scale]);
				if (nearestMode(modesByScale.get(scale), apex, toleranceDa / cfg.binWidth) >= 0) {
					persistent++;
				}
			}
			if (persistent < Math.min(cfg.kdeMinPersistentScales, bandwidths.length)) {
				continue;
			}

			int supportRadius = Math.max(1,
					(int) Math.ceil(Math.max(cfg.membershipFloorDa, 3.0 * bandwidths[0]) / cfg.binWidth));
			double targetSupport = sum(target, apex - supportRadius, apex + supportRadius);
			double decoySupport = sum(decoy, apex - supportRadius, apex + supportRadius);
			if (targetSupport - decoySupport < cfg.minPeakArea) {
				continue;
			}

			int reproducibleFiles = reproducibleFiles(apex, supportRadius, targetByFile, decoyByFile,
					bandwidths[0] / cfg.binWidth, cfg);
			if (!targetByFile.isEmpty() && reproducibleFiles < Math.min(cfg.kdeMinReproducibleFiles,
					targetByFile.size())) {
				continue;
			}
			double prominence = prominenceRatio(base, apex, supportRadius);
			if (cfg.kdeMinProminenceRatio > 0.0 && prominence < cfg.kdeMinProminenceRatio) {
				continue;
			}
			candidates.add(new Candidate(apex, persistent, reproducibleFiles, targetSupport, decoySupport,
					prominence));
		}

		candidates.sort((a, b) -> Double.compare(b.targetSupport - b.decoySupport,
				a.targetSupport - a.decoySupport));
		int evidenceQualified = candidates.size();

		// Resolution gate. Two modes closer together than a few standard
		// deviations are not two features: at this instrument's precursor
		// accuracy they cannot be told apart, so keeping both splits one peak
		// into several. Measured on the 12mix run, this is the dominant artefact
		// - five separate "peaks" were reported inside 0.021 Da where the fitted
		// sigma is 0.002 Da, and the same pattern repeated at +2.9 and +3.9 Da.
		// Such fragments pass every statistical test (reproducible in 9/9 files,
		// persistent across scales) because a systematic artefact is reproducible
		// by construction, so this has to be enforced geometrically.
		//
		// cfg.minPeakSeparationSigmas already existed for exactly this purpose but
		// only in the Gaussian-fitting branch, which the KDE path bypasses; it has
		// therefore never taken effect on a real run.
		//
		// Candidates arrive sorted by target excess, so the strongest mode of a
		// split group is the one kept. Absorbed modes are dropped rather than
		// merged: membership boundaries below are valleys between surviving
		// neighbours, so removing a neighbour widens the survivor's window over
		// exactly the range the dropped mode covered.
		int splitAbsorbed = 0;
		if (cfg.minPeakSeparationSigmas > 0.0 && candidates.size() > 1) {
			double[] centroidOf = new double[candidates.size()];
			double[] sigmaOf = new double[candidates.size()];
			for (int i = 0; i < candidates.size(); i++) {
				centroidOf[i] = refinedMode(base, candidates.get(i).apex, cfg);
				sigmaOf[i] = localWidth(base, candidates.get(i).apex, bandwidths[0], cfg);
			}
			List<Candidate> separated = new ArrayList<>(candidates.size());
			List<double[]> kept = new ArrayList<>(candidates.size());
			for (int i = 0; i < candidates.size(); i++) {
				boolean absorbed = false;
				for (double[] other : kept) {
					double needed = cfg.minPeakSeparationSigmas * Math.max(sigmaOf[i], other[1]);
					if (Math.abs(centroidOf[i] - other[0]) < needed) {
						absorbed = true;
						break;
					}
				}
				if (absorbed) {
					splitAbsorbed++;
					continue;
				}
				separated.add(candidates.get(i));
				kept.add(new double[] { centroidOf[i], sigmaOf[i] });
			}
			candidates = separated;
		}

		// Width gate, off by default. A mode much broader than the instrument's
		// own mass accuracy is an unresolved blend rather than a peak. The scale
		// is taken from the data - the median width of the surviving modes - so
		// this carries no assumption about the instrument.
		int tooBroad = 0;
		if (cfg.kdeMaxSigmaMultiple > 0.0 && candidates.size() > 2) {
			double[] widths = new double[candidates.size()];
			for (int i = 0; i < candidates.size(); i++) {
				widths[i] = localWidth(base, candidates.get(i).apex, bandwidths[0], cfg);
			}
			double[] sorted = widths.clone();
			Arrays.sort(sorted);
			double limit = cfg.kdeMaxSigmaMultiple * sorted[sorted.length / 2];
			List<Candidate> narrow = new ArrayList<>(candidates.size());
			for (int i = 0; i < candidates.size(); i++) {
				if (widths[i] > limit) {
					tooBroad++;
				} else {
					narrow.add(candidates.get(i));
				}
			}
			candidates = narrow;
		}

		int beforeCap = candidates.size();
		if (cfg.kdeMaxModes > 0 && candidates.size() > cfg.kdeMaxModes) {
			candidates = new ArrayList<>(candidates.subList(0, cfg.kdeMaxModes));
		}
		StringBuilder prominenceReport = new StringBuilder();
		if (!candidates.isEmpty()) {
			double[] proms = new double[candidates.size()];
			for (int i = 0; i < candidates.size(); i++) {
				proms[i] = candidates.get(i).prominenceRatio;
			}
			Arrays.sort(proms);
			prominenceReport.append("Delta-mass KDE: prominence ratio of retained modes, deciles:");
			for (int d = 0; d <= 10; d++) {
				int at = Math.min(proms.length - 1, (int) Math.round(d / 10.0 * (proms.length - 1)));
				prominenceReport.append(String.format(java.util.Locale.ROOT, " %d%%=%.2f", d * 10, proms[at]));
			}
		}

		String summary = "Delta-mass KDE: " + evidenceQualified
				+ " modes passed persistence/support/reproducibility; " + splitAbsorbed
				+ " absorbed as unresolved splits (<" + cfg.minPeakSeparationSigmas + " sigma apart)"
				+ (cfg.kdeMaxSigmaMultiple > 0.0 ? "; " + tooBroad + " dropped as over-broad" : "")
				+ "; " + candidates.size() + " retained"
				+ (beforeCap > candidates.size()
						? " after kdeMaxModes=" + cfg.kdeMaxModes + " TRUNCATED " + (beforeCap - candidates.size())
								+ " modes - raise kdeMaxModes, the cap is deciding the peak set"
						: " (cap " + cfg.kdeMaxModes + " not binding)");
		LOGGER.info(summary);
		if (progress != null) {
			progress.accept(summary);
			if (prominenceReport.length() > 0) {
				progress.accept(prominenceReport.toString());
			}
		}
		candidates.sort((a, b) -> Integer.compare(a.apex, b.apex));

		List<DeltaMassPeakFitter.Peak> peaks = new ArrayList<>(candidates.size());
		for (int i = 0; i < candidates.size(); i++) {
			Candidate candidate = candidates.get(i);
			int previous = i == 0 ? 0 : candidates.get(i - 1).apex;
			int next = i + 1 == candidates.size() ? bins - 1 : candidates.get(i + 1).apex;
			int cap = Math.max(1, (int) Math.ceil(cfg.kdeMaxMembershipDa / cfg.binWidth));
			int floor = Math.max(1, (int) Math.ceil(cfg.membershipFloorDa / cfg.binWidth));
			int left = valley(base, Math.max(previous + 1, candidate.apex - cap), candidate.apex, true);
			int right = valley(base, candidate.apex, Math.min(next - 1, candidate.apex + cap), false);
			left = Math.min(left, candidate.apex - floor);
			right = Math.max(right, candidate.apex + floor);
			left = Math.max(0, left);
			right = Math.min(bins - 1, right);

			double centroid = refinedMode(base, candidate.apex, cfg);
			double sigma = localWidth(base, candidate.apex, bandwidths[0], cfg);
			double amplitude = base[candidate.apex];
			double background = decoyBase == null ? 0.0 : decoyBase[candidate.apex];
			peaks.add(new DeltaMassPeakFitter.Peak(peaks.size(), centroid, sigma, amplitude, background,
					Double.NaN, 0, candidate.targetSupport, 0, -1,
					binCentre(left, cfg) - cfg.binWidth / 2.0,
					binCentre(right, cfg) + cfg.binWidth / 2.0,
					candidate.persistentScales, candidate.reproducibleFiles, "KDE"));
		}
		return peaks;
	}

	private static int reproducibleFiles(int apex, int supportRadius, Map<String, double[]> targets,
			Map<String, double[]> decoys, double bandwidthBins, DeltaMassPeakFitter.Config cfg) {
		int reproducible = 0;
		for (Map.Entry<String, double[]> entry : targets.entrySet()) {
			double[] target = entry.getValue();
			double[] decoy = decoys.getOrDefault(entry.getKey(), new double[target.length]);
			double targetSupport = sum(target, apex - supportRadius, apex + supportRadius);
			double decoySupport = sum(decoy, apex - supportRadius, apex + supportRadius);
			if (targetSupport < cfg.kdeMinFileSupport || targetSupport <= decoySupport) {
				continue;
			}
			double[] targetDensity = gaussianSmoothWindow(target, apex, bandwidthBins, supportRadius);
			double[] decoyDensity = gaussianSmoothWindow(decoy, apex, bandwidthBins, supportRadius);
			int centre = supportRadius;
			double centreExcess = targetDensity[centre] - decoyDensity[centre];
			double localMax = Double.NEGATIVE_INFINITY;
			for (int i = Math.max(0, centre - 2); i <= Math.min(targetDensity.length - 1, centre + 2); i++) {
				localMax = Math.max(localMax, targetDensity[i] - decoyDensity[i]);
			}
			if (centreExcess > 0.0 && centreExcess >= 0.8 * localMax) {
				reproducible++;
			}
		}
		return reproducible;
	}

	private static double[] gaussianSmoothWindow(double[] histogram, int apex, double bandwidthBins, int radius) {
		int from = apex - radius;
		int to = apex + radius;
		double[] out = new double[2 * radius + 1];
		int kernelRadius = Math.max(1, (int) Math.ceil(4.0 * bandwidthBins));
		double[] kernel = gaussianKernel(bandwidthBins, kernelRadius);
		for (int local = 0; local < out.length; local++) {
			int bin = from + local;
			double value = 0.0;
			for (int k = -kernelRadius; k <= kernelRadius; k++) {
				int source = bin - k;
				if (source >= 0 && source < histogram.length) {
					value += histogram[source] * kernel[k + kernelRadius];
				}
			}
			out[local] = value;
		}
		return out;
	}

	static double[] gaussianSmooth(double[] histogram, double bandwidthBins) {
		if (!(bandwidthBins > 0.0)) {
			return Arrays.copyOf(histogram, histogram.length);
		}
		int radius = Math.max(1, (int) Math.ceil(4.0 * bandwidthBins));
		double[] kernel = gaussianKernel(bandwidthBins, radius);
		double[] out = new double[histogram.length];
		for (int source = 0; source < histogram.length; source++) {
			if (histogram[source] == 0.0) {
				continue;
			}
			for (int k = -radius; k <= radius; k++) {
				int destination = source + k;
				if (destination >= 0 && destination < out.length) {
					out[destination] += histogram[source] * kernel[k + radius];
				}
			}
		}
		return out;
	}

	private static double[] gaussianKernel(double bandwidthBins, int radius) {
		double[] kernel = new double[2 * radius + 1];
		double total = 0.0;
		for (int k = -radius; k <= radius; k++) {
			double z = k / bandwidthBins;
			kernel[k + radius] = Math.exp(-0.5 * z * z);
			total += kernel[k + radius];
		}
		for (int i = 0; i < kernel.length; i++) {
			kernel[i] /= total;
		}
		return kernel;
	}

	private static int[] localMaxima(double[] density) {
		int[] temporary = new int[Math.max(1, density.length / 4)];
		int count = 0;
		for (int i = 1; i < density.length - 1; i++) {
			if (density[i] <= 0.0 || !(density[i] > density[i - 1] && density[i] >= density[i + 1])) {
				continue;
			}
			if (count == temporary.length) {
				temporary = Arrays.copyOf(temporary, temporary.length * 2);
			}
			temporary[count++] = i;
		}
		return Arrays.copyOf(temporary, count);
	}

	private static int nearestMode(int[] modes, int apex, double toleranceBins) {
		int position = Arrays.binarySearch(modes, apex);
		if (position >= 0) {
			return modes[position];
		}
		position = -position - 1;
		int best = -1;
		double distance = Double.POSITIVE_INFINITY;
		for (int candidate : new int[] { position - 1, position }) {
			if (candidate >= 0 && candidate < modes.length) {
				double current = Math.abs(modes[candidate] - apex);
				if (current < distance) {
					distance = current;
					best = modes[candidate];
				}
			}
		}
		return distance <= toleranceBins ? best : -1;
	}

	private static int valley(double[] density, int from, int to, boolean preferRightOnTie) {
		if (from > to) {
			return preferRightOnTie ? to : from;
		}
		int best = from;
		for (int i = from + 1; i <= to; i++) {
			if (density[i] < density[best]
					|| (preferRightOnTie && density[i] == density[best])) {
				best = i;
			}
		}
		return best;
	}

	private static double refinedMode(double[] density, int apex, DeltaMassPeakFitter.Config cfg) {
		double centre = binCentre(apex, cfg);
		if (apex <= 0 || apex >= density.length - 1) {
			return centre;
		}
		double left = density[apex - 1];
		double middle = density[apex];
		double right = density[apex + 1];
		double denominator = left - 2.0 * middle + right;
		if (Math.abs(denominator) < 1e-12) {
			return centre;
		}
		double fraction = 0.5 * (left - right) / denominator;
		fraction = Math.max(-0.5, Math.min(0.5, fraction));
		return centre + fraction * cfg.binWidth;
	}

	private static double localWidth(double[] density, int apex, double bandwidthDa,
			DeltaMassPeakFitter.Config cfg) {
		int radius = Math.max(1, (int) Math.ceil(3.0 * bandwidthDa / cfg.binWidth));
		double centre = binCentre(apex, cfg);
		double total = 0.0;
		double variance = 0.0;
		for (int bin = Math.max(0, apex - radius); bin <= Math.min(density.length - 1, apex + radius); bin++) {
			double weight = density[bin];
			double delta = binCentre(bin, cfg) - centre;
			total += weight;
			variance += weight * delta * delta;
		}
		return total > 0.0
				? Math.max(cfg.minSigmaBins * cfg.binWidth, Math.sqrt(variance / total))
				: Math.max(cfg.minSigmaBins * cfg.binWidth, bandwidthDa);
	}

	private static double sum(double[] values, int from, int to) {
		double total = 0.0;
		for (int i = Math.max(0, from); i <= Math.min(values.length - 1, to); i++) {
			total += values[i];
		}
		return total;
	}

	private static int binOf(double mass, DeltaMassPeakFitter.Config cfg) {
		return (int) Math.floor((mass - cfg.windowLower) / cfg.binWidth);
	}

	private static double binCentre(int bin, DeltaMassPeakFitter.Config cfg) {
		return cfg.windowLower + (bin + 0.5) * cfg.binWidth;
	}

	private static final class Candidate {
		final int apex;
		final int persistentScales;
		final int reproducibleFiles;
		final double targetSupport;
		final double decoySupport;
		/** Height above the higher flanking valley, as a fraction of apex height. */
		final double prominenceRatio;

		Candidate(int apex, int persistentScales, int reproducibleFiles, double targetSupport,
				double decoySupport, double prominenceRatio) {
			this.apex = apex;
			this.persistentScales = persistentScales;
			this.reproducibleFiles = reproducibleFiles;
			this.targetSupport = targetSupport;
			this.decoySupport = decoySupport;
			this.prominenceRatio = prominenceRatio;
		}
	}

	/**
	 * Topographic prominence of a mode, as a fraction of its own height: how far
	 * it rises above the higher of the two valleys flanking it.
	 * <p>
	 * A local maximum is not a peak. Between roughly +1.9 and +4 Da the 12mix
	 * delta-mass axis carries a broad elevated region, and mode-finding alone
	 * carved it into about twelve "peaks" per 0.06 Da - every wiggle on the hump
	 * became an entry, each passing the absolute support test (minPeakArea = 15)
	 * thousands of times over. Prominence is what distinguishes a peak standing on
	 * a background from a ripple in it, and it is the same measure PTM-Shepherd
	 * uses for its own peak picking.
	 *
	 * @return 1.0 for an isolated peak on a flat baseline, near 0 for a ripple.
	 */
	private static double prominenceRatio(double[] density, int apex, int searchRadius) {
		double height = density[apex];
		if (!(height > 0.0)) {
			return 0.0;
		}
		int lo = Math.max(0, apex - searchRadius);
		int hi = Math.min(density.length - 1, apex + searchRadius);
		double leftValley = height;
		for (int i = apex; i >= lo; i--) {
			leftValley = Math.min(leftValley, density[i]);
			if (density[i] > height) {
				break;
			}
		}
		double rightValley = height;
		for (int i = apex; i <= hi; i++) {
			rightValley = Math.min(rightValley, density[i]);
			if (density[i] > height) {
				break;
			}
		}
		return (height - Math.max(leftValley, rightValley)) / height;
	}
}
