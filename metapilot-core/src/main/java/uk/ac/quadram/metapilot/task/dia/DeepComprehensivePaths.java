package uk.ac.quadram.metapilot.task.dia;

import java.io.File;
import java.util.regex.Pattern;

/**
 * Folder layout for the Deep Comprehensive workflow (like DDA MAG per-run structure).
 * <p>
 * Result folder layout:
 * <ul>
 *   <li>{@code marker/} – Marker4 library search results (parquet, etc.)</li>
 *   <li>{@code focus/} – Focused library search results</li>
 *   <li>{@code training/} – Training dataset TSV and trained model</li>
 *   <li>{@code training/model/} – Trained model artifacts</li>
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

	/** Subfolder for training dataset and trained model. */
	public static final String FOLDER_TRAINING = "training";

	/** Subfolder under {@link #FOLDER_TRAINING} for the trained model. */
	public static final String FOLDER_MODEL = "model";

	/** Default training dataset filename under {@link #FOLDER_TRAINING}. */
	public static final String TRAINING_DATASET_TSV = "training_seed.tsv";

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
	 * Directory for training dataset and trained model.
	 */
	public static File getTrainingDir(File resultFolder) {
		return new File(resultFolder, FOLDER_TRAINING);
	}

	/**
	 * Directory for the trained model (under training).
	 */
	public static File getModelDir(File resultFolder) {
		return new File(getTrainingDir(resultFolder), FOLDER_MODEL);
	}

	/**
	 * Default file for the training dataset TSV.
	 */
	public static File getTrainingDatasetFile(File resultFolder) {
		return new File(getTrainingDir(resultFolder), TRAINING_DATASET_TSV);
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
	 * Ensure the training folder and model subfolder exist; creates parent dirs if needed.
	 *
	 * @return training dir
	 */
	public static File ensureTrainingFolders(File resultFolder) {
		File training = getTrainingDir(resultFolder);
		File model = getModelDir(resultFolder);
		if (!training.exists()) {
			training.mkdirs();
		}
		if (!model.exists()) {
			model.mkdirs();
		}
		return training;
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
