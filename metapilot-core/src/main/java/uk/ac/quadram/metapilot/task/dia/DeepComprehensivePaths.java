package uk.ac.quadram.metapilot.task.dia;

import java.io.File;
import java.util.regex.Pattern;

/**
 * Folder layout for the Seed Deep workflow (seed_deep result folder).
 * <p>
 * Result folder layout:
 * <ul>
 *   <li>{@code marker/} – Marker4 library search results (parquet, etc.)</li>
 *   <li>{@code focus/} – Focused library search results</li>
 *   <li>{@code data/} – Seed dataset TSV and optional model artifacts</li>
 *   <li>{@code data/model/} – Optional model artifacts (for future use)</li>
 *   <li>{@code <runId>/} – One folder per raw file: predict dataset and predictions for that run</li>
 * </ul>
 * Run IDs are derived from experiment name or raw file stem and sanitized for use as directory names.
 */
public final class DeepComprehensivePaths {

	private static final Pattern INVALID_PATH = Pattern.compile("[\\\\/:*?\"<>|]");

	private DeepComprehensivePaths() {}

	/** Subfolder for Marker4 search results. */
	public static final String FOLDER_MARKER = "marker";

	/** Subfolder for focused library search results. */
	public static final String FOLDER_FOCUS = "focus";

	/** Prefix for iterative search round subfolders under focus (round1, round2, ...). */
	public static final String FOLDER_ROUND_PREFIX = "round";

	/** Subfolder for seed dataset and optional model artifacts. */
	public static final String FOLDER_DATA = "data";

	/** Subfolder under {@link #FOLDER_DATA} for optional model artifacts. */
	public static final String FOLDER_MODEL = "model";

	/** Seed dataset filename under {@link #FOLDER_DATA}. */
	public static final String SEED_DATASET_TSV = "seed_dataset.tsv";

	/** Default predict dataset filename under each per-run folder. */
	public static final String PREDICT_DATASET_TSV = "predict_dataset.tsv";

	/** Default predictions filename under each per-run folder. */
	public static final String PREDICTIONS_TSV = "predictions.tsv";

	/**
	 * Directory for Marker4 search output.
	 */
	public static File getMarkerDir(File resultFolder) {
		return new File(resultFolder, FOLDER_MARKER);
	}

	/**
	 * Directory for focused library search output.
	 */
	public static File getFocusDir(File resultFolder) {
		return new File(resultFolder, FOLDER_FOCUS);
	}

	/**
	 * Directory for iterative search round {@code roundIndex} (1-based). Round 0 is focus/combined.
	 * E.g. round 1 = focus/round1, round 2 = focus/round2.
	 */
	public static File getRoundDir(File resultFolder, int roundIndex) {
		return new File(getFocusDir(resultFolder), FOLDER_ROUND_PREFIX + roundIndex);
	}

	/**
	 * Directory for seed dataset and optional model (data folder).
	 */
	public static File getDataDir(File resultFolder) {
		return new File(resultFolder, FOLDER_DATA);
	}

	/**
	 * Directory for the optional model artifacts (under data).
	 */
	public static File getModelDir(File resultFolder) {
		return new File(getDataDir(resultFolder), FOLDER_MODEL);
	}

	/**
	 * Default file for the seed dataset TSV.
	 */
	public static File getSeedDatasetFile(File resultFolder) {
		return new File(getDataDir(resultFolder), SEED_DATASET_TSV);
	}

	/**
	 * Per-raw-file directory for predict dataset and predictions.
	 * Uses a sanitized run ID (experiment name or file stem).
	 *
	 * @param resultFolder result root
	 * @param runId        run identifier (e.g. experiment name or raw file name without path/extension)
	 * @return directory for this run (e.g. resultFolder/Exp1/)
	 */
	public static File getPerRunDir(File resultFolder, String runId) {
		return new File(resultFolder, sanitizeRunId(runId));
	}

	/**
	 * Predict dataset TSV for a given run.
	 */
	public static File getPredictDatasetFile(File resultFolder, String runId) {
		return new File(getPerRunDir(resultFolder, runId), PREDICT_DATASET_TSV);
	}

	/**
	 * Predictions TSV for a given run.
	 */
	public static File getPredictionsFile(File resultFolder, String runId) {
		return new File(getPerRunDir(resultFolder, runId), PREDICTIONS_TSV);
	}

	/**
	 * Sanitize a string for use as a directory name (replace invalid path chars with underscore).
	 */
	public static String sanitizeRunId(String runId) {
		if (runId == null || runId.isEmpty()) {
			return "run";
		}
		String s = INVALID_PATH.matcher(runId.trim()).replaceAll("_");
		return s.isEmpty() ? "run" : s;
	}

	/**
	 * Derive run ID from raw file path (filename without extension). Use when experiment name is not available.
	 */
	public static String runIdFromRawFilePath(String rawFilePath) {
		if (rawFilePath == null || rawFilePath.isEmpty()) {
			return "run";
		}
		int lastSlash = Math.max(rawFilePath.lastIndexOf('/'), rawFilePath.lastIndexOf('\\'));
		String name = lastSlash >= 0 ? rawFilePath.substring(lastSlash + 1) : rawFilePath;
		int dot = name.lastIndexOf('.');
		String stem = dot > 0 ? name.substring(0, dot) : name;
		return sanitizeRunId(stem);
	}

	/**
	 * Ensure the data folder and optional model subfolder exist; creates parent dirs if needed.
	 *
	 * @return data dir
	 */
	public static File ensureDataFolders(File resultFolder) {
		File data = getDataDir(resultFolder);
		File model = getModelDir(resultFolder);
		if (!data.exists()) {
			data.mkdirs();
		}
		if (!model.exists()) {
			model.mkdirs();
		}
		return data;
	}

	/**
	 * Ensure the per-run folder exists for the given run ID.
	 *
	 * @return per-run dir
	 */
	public static File ensurePerRunFolder(File resultFolder, String runId) {
		File dir = getPerRunDir(resultFolder, runId);
		if (!dir.exists()) {
			dir.mkdirs();
		}
		return dir;
	}
}
