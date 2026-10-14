package uk.ac.quadram.metapilot.core.opensearch;

import static org.junit.jupiter.api.Assertions.assertNotNull;
import static org.junit.jupiter.api.Assertions.assertNull;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.util.ArrayList;
import java.util.List;

import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Test;

class DeltaMassKdePeakFinderTest {

	private static DeltaMassPeakFitter.Config config() {
		DeltaMassPeakFitter.Config cfg = new DeltaMassPeakFitter.Config();
		cfg.windowLower = 0.5;
		cfg.windowUpper = 3.5;
		cfg.binWidth = 0.001;
		cfg.kdeBandwidthDa = 0.002;
		cfg.kdeScaleMultipliers = "1/2/4";
		cfg.kdeMinPersistentScales = 2;
		cfg.kdeMinReproducibleFiles = 2;
		cfg.kdeMinFileSupport = 3.0;
		cfg.minPeakArea = 8;
		return cfg;
	}

	private static DeltaMassPeakFitter.Peak nearest(DeltaMassPeakFitter.Result result, double mass) {
		DeltaMassPeakFitter.Peak best = null;
		double distance = Double.POSITIVE_INFINITY;
		for (DeltaMassPeakFitter.Peak peak : result.peaks) {
			double current = Math.abs(peak.centroid - mass);
			if (current < distance) {
				distance = current;
				best = peak;
			}
		}
		return distance <= 0.02 ? best : null;
	}

	@Test
	@DisplayName("KDE keeps a mode reproduced across files and rejects a one-file artefact")
	void fileReproducibilityIsPartOfDiscovery() {
		List<Double> masses = new ArrayList<>();
		List<String> files = new ArrayList<>();
		for (String file : new String[] { "a", "b", "c" }) {
			for (int i = 0; i < 40; i++) {
				masses.add(1.0 + ((i % 13) - 6) * 0.0007);
				files.add(file);
			}
		}
		for (int i = 0; i < 60; i++) {
			masses.add(3.0 + ((i % 9) - 4) * 0.0005);
			files.add("a");
		}

		DeltaMassPeakFitter.Result result = fit(masses, files, allTargets(masses.size()), config());
		DeltaMassPeakFitter.Peak shared = nearest(result, 1.0);
		assertNotNull(shared);
		assertTrue(shared.persistentScales >= 2);
		assertTrue(shared.reproducibleFiles >= 2);
		assertNull(nearest(result, 3.0), "a strong mode confined to one file is not reproducible");
	}

	@Test
	@DisplayName("Target-excess KDE suppresses a balanced target-decoy mode")
	void balancedTargetDecoyModeIsNotAModificationHypothesis() {
		List<Double> masses = new ArrayList<>();
		List<String> files = new ArrayList<>();
		List<Boolean> labels = new ArrayList<>();
		for (String file : new String[] { "a", "b", "c" }) {
			for (int i = 0; i < 30; i++) {
				double jitter = ((i % 11) - 5) * 0.0006;
				masses.add(1.0 + jitter);
				files.add(file);
				labels.add(true);
			}
			for (int i = 0; i < 30; i++) {
				double jitter = ((i % 11) - 5) * 0.0006;
				masses.add(2.0 + jitter);
				files.add(file);
				labels.add(true);
				masses.add(2.0 + jitter);
				files.add(file);
				labels.add(false);
			}
		}
		boolean[] targets = new boolean[labels.size()];
		for (int i = 0; i < targets.length; i++) {
			targets[i] = labels.get(i);
		}

		DeltaMassPeakFitter.Result result = fit(masses, files, targets, config());
		assertNotNull(nearest(result, 1.0));
		assertNull(nearest(result, 2.0), "equal target and decoy density has no target excess");
	}

	private static DeltaMassPeakFitter.Result fit(List<Double> massList, List<String> fileList,
			boolean[] targets, DeltaMassPeakFitter.Config cfg) {
		double[] masses = massList.stream().mapToDouble(Double::doubleValue).toArray();
		String[] files = fileList.toArray(new String[0]);
		boolean[] discoveryMask = targets.clone();
		double[] weights = new double[masses.length];
		java.util.Arrays.fill(weights, 1.0);
		return DeltaMassPeakFitter.fit(masses, discoveryMask, weights, targets, files, cfg, null);
	}

	private static boolean[] allTargets(int length) {
		boolean[] targets = new boolean[length];
		java.util.Arrays.fill(targets, true);
		return targets;
	}
}
