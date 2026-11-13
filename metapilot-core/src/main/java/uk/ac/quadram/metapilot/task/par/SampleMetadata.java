package uk.ac.quadram.metapilot.task.par;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Locale;
import java.util.Map;

import org.json.JSONArray;
import org.json.JSONObject;

/**
 * One metadata row per acquisition, shared by every workflow parameter file
 * and by the report's {@code metainfo.tsv}.
 * <p>
 * Row order is {@code sample name, label type, age, gender, condition,
 * treatment, time point, custom 1, custom 2, custom 3}. The sample name is
 * the experiment name. Condition is what the HTML report groups on.
 */
public final class SampleMetadata {

	public static final String JSON_KEY = "metadata";

	private static final String[] FIELD_TITLES = {
			"Label Type", "Age", "Gender", "Condition", "Treatment", "Time Point",
			"Custom 1", "Custom 2", "Custom 3"
	};
	private static final String[] FIELD_KEYS = {
			"label_type", "age", "gender", "condition", "treatment", "time_point",
			"custom_field_1", "custom_field_2", "custom_field_3"
	};

	private SampleMetadata() {
	}

	public static JSONArray toJson(List<String[]> rows) {
		JSONArray metadata = new JSONArray();
		if (rows == null) {
			return metadata;
		}
		for (String[] row : rows) {
			if (row != null) {
				metadata.put(rowToJson(row));
			}
		}
		return metadata;
	}

	/**
	 * Writes the metadata already stored on {@code metadata}. Column titles choose
	 * the JSON field, so a visible Condition column is saved as {@code condition}
	 * rather than as whichever position it occupied in the table.
	 */
	public static JSONArray toJson(MetaData metadata) {
		JSONArray array = new JSONArray();
		if (metadata == null || metadata.getMetaInfo() == null) {
			return array;
		}
		String[][] info = metadata.getMetaInfo();
		String[] titles = metadata.getMetaTypeTitles();
		boolean titled = titles != null && titles.length > 0;
		for (int i = 0; i < info.length; i++) {
			String[] stored = blankRow();
			stored[0] = slotName(metadata, i);
			String[] values = info[i] == null ? new String[0] : info[i];
			if (titled) {
				for (int c = 0; c < values.length && c < titles.length; c++) {
					int field = fieldIndex(titles[c]);
					if (field >= 0) {
						stored[field] = values[c] == null ? "" : values[c];
					}
				}
			} else {
				// Files written before column titles were stored used this order.
				if (values.length > 0) {
					stored[2] = values[0] == null ? "" : values[0];
				}
				if (values.length > 1) {
					stored[3] = values[1] == null ? "" : values[1];
				}
				if (values.length > 2) {
					stored[4] = values[2] == null ? "" : values[2];
				}
			}
			array.put(rowToJson(stored));
		}
		return array;
	}

	/**
	 * @return the metadata rows, an empty list when the array is present but empty,
	 *         or null when the file has no metadata section
	 */
	public static List<String[]> fromJson(JSONObject json) {
		if (json == null || !json.has(JSON_KEY)) {
			return null;
		}
		JSONArray array = json.optJSONArray(JSON_KEY);
		if (array == null) {
			return new ArrayList<>();
		}
		boolean standard = false;
		for (int i = 0; i < array.length(); i++) {
			JSONObject row = array.optJSONObject(i);
			if (row != null && (row.has("sample_name") || row.has("label_type") || row.has("time_point")
					|| row.has("custom_field_1"))) {
				standard = true;
				break;
			}
		}
		if (standard) {
			List<String[]> rows = new ArrayList<>();
			for (int i = 0; i < array.length(); i++) {
				JSONObject row = array.optJSONObject(i);
				if (row != null) {
					rows.add(readStandardRow(row));
				}
			}
			return rows;
		}
		return readPositional(array);
	}

	/**
	 * Copies metadata rows onto {@code metadata} so {@link MetaData#exportMetadata}
	 * can write them. Rows line up with acquisitions when the counts match, which
	 * is how the metadata table is kept in step with the file list. Otherwise a
	 * row is used when its sample name is the experiment name.
	 */
	public static void apply(MetaData metadata, List<String[]> rows) {
		if (metadata == null || rows == null || rows.isEmpty()) {
			return;
		}
		List<Integer> active = new ArrayList<>();
		for (int field = 1; field <= FIELD_TITLES.length; field++) {
			if (columnHasValue(rows, field)) {
				active.add(field);
			}
		}
		if (active.isEmpty()) {
			return;
		}
		int slots = slotCount(metadata);
		if (slots <= 0) {
			slots = rows.size();
		}
		String[] titles = new String[active.size()];
		for (int c = 0; c < active.size(); c++) {
			titles[c] = FIELD_TITLES[active.get(c) - 1];
		}
		String[][] info = new String[slots][active.size()];
		Map<String, String[]> byName = new LinkedHashMap<>();
		for (String[] row : rows) {
			if (row != null && row.length > 0 && row[0] != null && !row[0].isBlank()) {
				byName.putIfAbsent(row[0], row);
			}
		}
		boolean byIndex = rows.size() == slots;
		for (int i = 0; i < slots; i++) {
			String[] row = null;
			if (byIndex) {
				row = rows.get(i);
			} else {
				String name = slotName(metadata, i);
				if (name != null && !name.isBlank()) {
					row = byName.get(name);
				}
			}
			for (int c = 0; c < active.size(); c++) {
				int field = active.get(c);
				String value = "";
				if (row != null && field < row.length && row[field] != null) {
					value = row[field];
				}
				info[i][c] = value;
			}
		}
		metadata.setMetaTypeCount(active.size());
		metadata.setMetaTypeTitles(titles);
		metadata.setMetaInfo(info);
	}

	/**
	 * Writes {@code resultDir/metainfo.tsv} when the acquisition metadata has
	 * columns to report. Returns null when there is nothing to group on.
	 */
	public static File reportFile(File resultDir, MetaData metadata) throws IOException {
		if (resultDir == null || metadata == null || metadata.getMetaTypeCount() <= 0) {
			return null;
		}
		if (!resultDir.exists() && !resultDir.mkdirs()) {
			return null;
		}
		File file = new File(resultDir, "metainfo.tsv");
		metadata.exportMetadata(file);
		return file;
	}

	private static JSONObject rowToJson(String[] row) {
		JSONObject meta = new JSONObject();
		meta.put("sample_name", value(row, 0));
		for (int i = 0; i < FIELD_KEYS.length; i++) {
			meta.put(FIELD_KEYS[i], value(row, i + 1));
		}
		return meta;
	}

	private static String[] readStandardRow(JSONObject row) {
		String[] metadataRow = blankRow();
		metadataRow[0] = row.optString("sample_name", "");
		for (int i = 0; i < FIELD_KEYS.length; i++) {
			metadataRow[i + 1] = row.optString(FIELD_KEYS[i], "");
		}
		return metadataRow;
	}

	private static List<String[]> readPositional(JSONArray array) {
		int slots = array.length();
		for (int i = 0; i < array.length(); i++) {
			JSONObject row = array.optJSONObject(i);
			if (row != null && row.has("sample_index")) {
				slots = Math.max(slots, row.optInt("sample_index", i) + 1);
			}
		}
		List<String[]> rows = new ArrayList<>();
		for (int i = 0; i < slots; i++) {
			rows.add(blankRow());
		}
		for (int i = 0; i < array.length(); i++) {
			JSONObject row = array.optJSONObject(i);
			if (row == null) {
				continue;
			}
			int index = row.has("sample_index") ? row.optInt("sample_index", i) : i;
			if (index < 0 || index >= rows.size()) {
				continue;
			}
			String[] metadataRow = rows.get(index);
			metadataRow[2] = row.optString("age", "");
			metadataRow[3] = row.optString("gender", "");
			metadataRow[4] = row.optString("condition", "");
		}
		return rows;
	}

	private static boolean columnHasValue(List<String[]> rows, int field) {
		for (String[] row : rows) {
			if (row != null && field < row.length && row[field] != null && !row[field].isBlank()) {
				return true;
			}
		}
		return false;
	}

	private static int fieldIndex(String title) {
		if (title == null || title.isBlank()) {
			return -1;
		}
		String key = title.trim().toLowerCase(Locale.ROOT);
		if ("bmi".equals(key) || "custom 1".equals(key)) {
			return 7;
		}
		for (int i = 0; i < FIELD_TITLES.length; i++) {
			if (FIELD_TITLES[i].equalsIgnoreCase(title.trim())) {
				return i + 1;
			}
		}
		return -1;
	}

	private static int slotCount(MetaData metadata) {
		String[] experiments = metadata.getExpNames();
		int experimentsCount = experiments == null ? 0 : experiments.length;
		String[] labels = metadata.getLabelTitle();
		if (labels != null && labels.length > 0 && experimentsCount > 0) {
			return experimentsCount * labels.length;
		}
		if (experimentsCount > 0) {
			return experimentsCount;
		}
		String[] rawFiles = metadata.getRawFiles();
		return rawFiles == null ? 0 : rawFiles.length;
	}

	private static String slotName(MetaData metadata, int index) {
		String[] labelNames = metadata.getLabelExpNames();
		if (labelNames != null && index >= 0 && index < labelNames.length && labelNames[index] != null
				&& !labelNames[index].isBlank()) {
			return labelNames[index];
		}
		String[] experiments = metadata.getExpNames();
		String[] labels = metadata.getLabelTitle();
		if (labels != null && labels.length > 0 && experiments != null && experiments.length > 0) {
			int experimentIndex = index / labels.length;
			int labelIndex = index % labels.length;
			if (experimentIndex < experiments.length && experiments[experimentIndex] != null) {
				return experiments[experimentIndex] + "_" + labels[labelIndex];
			}
		}
		if (experiments != null && index >= 0 && index < experiments.length && experiments[index] != null) {
			return experiments[index];
		}
		return "";
	}

	private static String[] blankRow() {
		return new String[] { "", "", "", "", "", "", "", "", "", "" };
	}

	private static String value(String[] row, int index) {
		return row != null && index < row.length && row[index] != null ? row[index] : "";
	}
}
