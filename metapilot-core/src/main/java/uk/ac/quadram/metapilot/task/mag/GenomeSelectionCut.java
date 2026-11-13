package uk.ac.quadram.metapilot.task.mag;

import java.util.Locale;

import org.json.JSONObject;

/**
 * PSM-coverage cut adopted from MetaLab-MAG (J Proteome Res 2023).
 *
 * The cut runs on MetaPilot's existing ranking. It does not rank by PSM cover,
 * and it is not the MetaLab-MAG workflow.
 */
public final class GenomeSelectionCut {

	public static final String STRATEGY_FOCUS = "Focus";
	public static final String STRATEGY_COMPREHENSIVE = "Comprehensive";
	public static final String STRATEGY_PSM_COVERAGE = "PSM coverage";

	public static final String STOP_COVERAGE = "COVERAGE";
	public static final String STOP_SEEDS_ABOVE_THRESHOLD = "SEEDS_ABOVE_THRESHOLD";
	public static final String STOP_STEP_GUARD = "STEP_GUARD";
	public static final String STOP_GENOME_CAP = "GENOME_CAP";
	public static final String STOP_END_OF_RANKING = "END_OF_RANKING";
	public static final String STOP_KNEE = "KNEE";
	public static final String STOP_FALLBACK_MIN1 = "FALLBACK_MIN1";
	public static final String STOP_FULL_COVERAGE = "FULL_COVERAGE";

	public static final double DEFAULT_PSM_COVERAGE = 0.95;
	/** Fraction of marker PSMs. The GUI shows this as 0.1%. */
	public static final double DEFAULT_STEP_GUARD = 0.001;
	/** Inclusive upper bound: 0.2% of marker PSMs. */
	public static final double MAX_STEP_GUARD = 0.002;
	/** Smallest GUI step: 0.01 percentage points, stored as this fraction. */
	public static final double STEP_GUARD_QUANTUM = 0.0001;
	public static final int DEFAULT_MAX_GENOMES = 300;

	public static final String KEY_PSM_COVERAGE = "genome_selection_psm_coverage";
	public static final String KEY_STEP_GUARD = "genome_selection_step_guard";
	public static final String KEY_MAX_GENOMES = "genome_selection_max_genomes";

	private GenomeSelectionCut() {
	}

	public static final class Decision {
		public final int selectedCount;
		public final double coverageReached;
		public final String stopReason;

		public Decision(int selectedCount, double coverageReached, String stopReason) {
			this.selectedCount = selectedCount;
			this.coverageReached = coverageReached;
			this.stopReason = stopReason;
		}
	}

	/**
	 * Walks the in-ranking genomes (seeds, then greedy order). Indexes are 0-based.
	 * Genomes listed after the greedy cap are not part of these arrays.
	 */
	public static Decision decide(double[] coverage, int[] increment, int totalPsms, int seedCount, double threshold,
			double stepGuard, int maxGenomes) {
		int n = coverage == null ? 0 : coverage.length;
		if (n == 0) {
			return new Decision(0, Double.NaN, STOP_END_OF_RANKING);
		}
		int guard = n - 1;
		if (stepGuard > 0) {
			int lastOverOne = -1;
			int lastOverStep = -1;
			double total = totalPsms;
			for (int i = 0; i < n; i++) {
				int gained = increment == null || i >= increment.length ? 0 : increment[i];
				if (gained > 1) {
					lastOverOne = i;
				}
				if (total > 0 && gained / total > stepGuard) {
					lastOverStep = i;
				}
			}
			guard = Math.min(lastOverOne < 0 ? 0 : lastOverOne, lastOverStep < 0 ? 0 : lastOverStep);
		}
		int minCheck = Math.max(seedCount, 1);
		for (int i = 0; i < n; i++) {
			int count = i + 1;
			if (count < minCheck) {
				continue;
			}
			if (coverage[i] >= threshold) {
				String why = count == seedCount ? STOP_SEEDS_ABOVE_THRESHOLD : STOP_COVERAGE;
				return new Decision(count, coverage[i], why);
			}
			if (stepGuard > 0 && i >= guard) {
				return new Decision(count, coverage[i], STOP_STEP_GUARD);
			}
			if (count >= maxGenomes) {
				return new Decision(count, coverage[i], STOP_GENOME_CAP);
			}
		}
		return new Decision(n, coverage[n - 1], STOP_END_OF_RANKING);
	}

	public static boolean isPsmCoverage(String strategy) {
		return strategy != null && STRATEGY_PSM_COVERAGE.equalsIgnoreCase(strategy.trim());
	}

	public static boolean isComprehensive(String strategy) {
		return strategy != null && STRATEGY_COMPREHENSIVE.equalsIgnoreCase(strategy.trim());
	}

	/**
	 * Blank becomes Focus. Anything else must be Focus, Comprehensive, or PSM coverage.
	 */
	public static String requireStrategy(String strategy) {
		if (strategy == null || strategy.trim().isEmpty()) {
			return STRATEGY_FOCUS;
		}
		String trimmed = strategy.trim();
		if (STRATEGY_FOCUS.equalsIgnoreCase(trimmed)) {
			return STRATEGY_FOCUS;
		}
		if (STRATEGY_COMPREHENSIVE.equalsIgnoreCase(trimmed)) {
			return STRATEGY_COMPREHENSIVE;
		}
		if (STRATEGY_PSM_COVERAGE.equalsIgnoreCase(trimmed)) {
			return STRATEGY_PSM_COVERAGE;
		}
		throw new IllegalArgumentException("genome_selection_strategy must be Focus, Comprehensive, or \""
				+ STRATEGY_PSM_COVERAGE + "\" (got " + strategy + ")");
	}

	public static void requireCoverage(double coverage) {
		if (!(coverage > 0.0 && coverage <= 1.0)) {
			throw new IllegalArgumentException(
					"genome_selection_psm_coverage must be > 0 and <= 1 (got " + coverage + ")");
		}
	}

	public static void requireStepGuard(double stepGuard) {
		if (!(stepGuard >= 0.0 && stepGuard <= MAX_STEP_GUARD + 1e-12)) {
			throw new IllegalArgumentException(
					"genome_selection_step_guard must be from 0 to 0.002 (0 to 0.2%; got " + stepGuard + ")");
		}
	}

	/** Spinner position for a stored fraction. 0.001 is step 10, shown as 0.1%. */
	public static int stepGuardPercentSteps(double fraction) {
		requireStepGuard(fraction);
		int steps = (int) Math.round(fraction / STEP_GUARD_QUANTUM);
		int max = (int) Math.round(MAX_STEP_GUARD / STEP_GUARD_QUANTUM);
		if (steps < 0 || steps > max) {
			throw new IllegalArgumentException(
					"genome_selection_step_guard must be a multiple of 0.0001 from 0 to 0.002 (got " + fraction + ")");
		}
		return steps;
	}

	/** Stored fraction for a spinner position. Step 10 is 0.001. */
	public static double stepGuardFromPercentSteps(int steps) {
		int max = (int) Math.round(MAX_STEP_GUARD / STEP_GUARD_QUANTUM);
		if (steps < 0 || steps > max) {
			throw new IllegalArgumentException(
					"step guard percentage must be from 0 to 0.2 in steps of 0.01 (got " + (steps / 100.0) + "%)");
		}
		return steps / 10000.0;
	}

	/** Parses a percentage such as {@code 0.1} or {@code 0.1%} into a spinner position. */
	public static int parseStepGuardPercent(String text) {
		if (text == null) {
			throw new IllegalArgumentException("step guard percentage is empty");
		}
		String trimmed = text.trim();
		if (trimmed.endsWith("%")) {
			trimmed = trimmed.substring(0, trimmed.length() - 1).trim();
		}
		if (trimmed.isEmpty()) {
			throw new IllegalArgumentException("step guard percentage is empty");
		}
		double percent;
		try {
			percent = Double.parseDouble(trimmed);
		} catch (NumberFormatException e) {
			throw new IllegalArgumentException("step guard percentage must be a number (got " + text + ")");
		}
		return stepGuardPercentSteps(percent / 100.0);
	}

	/** {@code 10} is {@code 0.1%}; {@code 1} is {@code 0.01%}. */
	public static String formatStepGuardPercent(int steps) {
		stepGuardFromPercentSteps(steps);
		int hundredths = steps;
		if (hundredths % 10 == 0) {
			return (hundredths / 100) + "." + (hundredths / 10 % 10) + "%";
		}
		return String.format(Locale.ROOT, "%d.%02d%%", hundredths / 100, hundredths % 100);
	}

	public static void requireMaxGenomes(int maxGenomes) {
		if (maxGenomes < 1) {
			throw new IllegalArgumentException(
					"genome_selection_max_genomes must be >= 1 (got " + maxGenomes + ")");
		}
	}

	/**
	 * Non-default coverage settings are stored with every parameter file so the
	 * controls can be restored, and they apply only when the strategy is PSM coverage.
	 */
	public static String unusedCoverageKeysWarning(String strategy, double coverage, double stepGuard, int maxGenomes) {
		if (isPsmCoverage(strategy)) {
			return null;
		}
		boolean nonDefault = Math.abs(coverage - DEFAULT_PSM_COVERAGE) > 1e-12
				|| Math.abs(stepGuard - DEFAULT_STEP_GUARD) > 1e-12
				|| maxGenomes != DEFAULT_MAX_GENOMES;
		if (!nonDefault) {
			return null;
		}
		String shown = strategy == null || strategy.trim().isEmpty() ? STRATEGY_FOCUS : strategy.trim();
		return "genome_selection_psm_coverage, genome_selection_step_guard and genome_selection_max_genomes "
				+ "are used only when genome_selection_strategy is \"" + STRATEGY_PSM_COVERAGE
				+ "\". The current strategy is " + shown + ", so those values are saved and ignored.";
	}

	public static String missedTargetMessage(double target, double reached, String stopReason, int maxGenomes) {
		if (Double.isNaN(target) || Double.isNaN(reached) || reached + 1e-12 >= target) {
			return "";
		}
		if (STOP_STEP_GUARD.equals(stopReason)) {
			return String.format(Locale.ROOT,
					"Reached %.1f%% of the %.0f%% target: stopped by the step guard. Lower or disable it to continue.",
					reached * 100.0, target * 100.0);
		}
		if (STOP_GENOME_CAP.equals(stopReason)) {
			return String.format(Locale.ROOT, "Stopped at Max genomes (%d) with %.1f%% coverage.", maxGenomes,
					reached * 100.0);
		}
		if (STOP_END_OF_RANKING.equals(stopReason)) {
			return String.format(Locale.ROOT, "Reached %.1f%% of the %.0f%% target: the ranking ended.",
					reached * 100.0, target * 100.0);
		}
		return "";
	}

	public static void writeAdditionalKeys(JSONObject additional, double coverage, double stepGuard, int maxGenomes) {
		additional.put(KEY_PSM_COVERAGE, coverage);
		additional.put(KEY_STEP_GUARD, stepGuard);
		additional.put(KEY_MAX_GENOMES, maxGenomes);
	}

	/** Unknown strategy text becomes Focus. Blank still becomes Focus, without a warning. */
	public static String strategyOrDefault(String strategy) {
		try {
			return requireStrategy(strategy);
		} catch (IllegalArgumentException e) {
			return STRATEGY_FOCUS;
		}
	}

	/** Null when {@code strategy} is blank or one of the three known names. */
	public static String strategyReplacementWarning(String strategy) {
		try {
			requireStrategy(strategy);
			return null;
		} catch (IllegalArgumentException e) {
			return replacementWarning("genome_selection_strategy", strategy, STRATEGY_FOCUS);
		}
	}

	/**
	 * Reads the three optional keys. Missing keys keep the supplied defaults.
	 * A present value that is not a number or is outside its range is replaced
	 * with that default. {@link LoadedCoverage#replacementWarning} lists those
	 * replacements.
	 */
	public static LoadedCoverage readAdditionalKeys(JSONObject additional, double coverage, double stepGuard,
			int maxGenomes) {
		boolean any = false;
		StringBuilder warnings = new StringBuilder();
		if (additional != null && additional.has(KEY_PSM_COVERAGE)) {
			any = true;
			coverage = acceptedDouble(additional, KEY_PSM_COVERAGE, coverage, GenomeSelectionCut::requireCoverage,
					warnings);
		}
		if (additional != null && additional.has(KEY_STEP_GUARD)) {
			any = true;
			stepGuard = acceptedDouble(additional, KEY_STEP_GUARD, stepGuard, GenomeSelectionCut::requireStepGuard,
					warnings);
		}
		if (additional != null && additional.has(KEY_MAX_GENOMES)) {
			any = true;
			maxGenomes = acceptedInt(additional, KEY_MAX_GENOMES, maxGenomes, warnings);
		}
		String warning = warnings.length() == 0 ? null : warnings.toString().trim();
		return new LoadedCoverage(any, coverage, stepGuard, maxGenomes, warning);
	}

	private static double acceptedDouble(JSONObject additional, String key, double fallback,
			java.util.function.DoubleConsumer check, StringBuilder warnings) {
		double stored;
		try {
			stored = requiredDouble(additional, key);
			check.accept(stored);
			return stored;
		} catch (RuntimeException e) {
			Object got = additional.has(key) ? additional.get(key) : storedText(e);
			appendWarning(warnings, replacementWarning(key, got, fallback));
			return fallback;
		}
	}

	private static int acceptedInt(JSONObject additional, String key, int fallback, StringBuilder warnings) {
		try {
			int stored = requiredInt(additional, key);
			requireMaxGenomes(stored);
			return stored;
		} catch (RuntimeException e) {
			Object got = additional.has(key) ? additional.get(key) : storedText(e);
			appendWarning(warnings, replacementWarning(key, got, fallback));
			return fallback;
		}
	}

	private static String storedText(RuntimeException e) {
		String message = e.getMessage();
		return message == null ? "" : message;
	}

	private static void appendWarning(StringBuilder warnings, String warning) {
		if (warnings.length() > 0) {
			warnings.append(' ');
		}
		warnings.append(warning);
	}

	private static String replacementWarning(String key, Object got, Object used) {
		return key + " " + got + " is not a usable value and was replaced with " + used + ".";
	}

	public static final class LoadedCoverage {
		public final boolean anyKey;
		public final double coverage;
		public final double stepGuard;
		public final int maxGenomes;
		/** Set when a stored value was replaced with its default. */
		public final String replacementWarning;

		public LoadedCoverage(boolean anyKey, double coverage, double stepGuard, int maxGenomes,
				String replacementWarning) {
			this.anyKey = anyKey;
			this.coverage = coverage;
			this.stepGuard = stepGuard;
			this.maxGenomes = maxGenomes;
			this.replacementWarning = replacementWarning;
		}
	}

	private static double requiredDouble(JSONObject object, String key) {
		Object value = object.get(key);
		try {
			double parsed = Double.parseDouble(String.valueOf(value));
			if (Double.isNaN(parsed) || Double.isInfinite(parsed)) {
				throw new NumberFormatException(String.valueOf(value));
			}
			return parsed;
		} catch (NumberFormatException e) {
			throw new IllegalArgumentException(key + " must be a number (got " + value + ")");
		}
	}

	private static int requiredInt(JSONObject object, String key) {
		Object value = object.get(key);
		try {
			double parsed = Double.parseDouble(String.valueOf(value));
			if (parsed != Math.rint(parsed) || parsed > Integer.MAX_VALUE || parsed < Integer.MIN_VALUE) {
				throw new NumberFormatException(String.valueOf(value));
			}
			return (int) parsed;
		} catch (NumberFormatException e) {
			throw new IllegalArgumentException(key + " must be an integer (got " + value + ")");
		}
	}

	public static String highTargetWarning(double target, double stepGuard) {
		if (target + 1e-12 < 0.97) {
			return "";
		}
		if (stepGuard > 0) {
			return "The step guard is likely to stop selection before this target.";
		}
		return "This can add hundreds of genomes per file on complex samples, and may hit Max genomes.";
	}
}
