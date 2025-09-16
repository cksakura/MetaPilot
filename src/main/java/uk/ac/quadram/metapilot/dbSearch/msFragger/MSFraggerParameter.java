/**
 * 
 */
package uk.ac.quadram.metapilot.dbSearch.msFragger;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Arrays;
import java.util.HashMap;

import uk.ac.quadram.metapilot.core.enzyme.Enzyme;

/**
 * @author Kai Cheng
 *
 */
public class MSFraggerParameter {

	private static final File params = new File("Resources\\MSFragger\\fragger.params");

	private HashMap<String, MSFraggerParItem> valueMap;
	private HashMap<String, String> variModMap;

	private static final int high_high = 0;
	private static final int high_low = 1;

	private MSFraggerParameter(HashMap<String, MSFraggerParItem> valueMap) {
		this.valueMap = valueMap;
	}

	public static MSFraggerParameter parse() throws IOException {
		if (params.exists()) {
			return parse(params);
		} else {
			return null;
		}
	}

	public static MSFraggerParameter parse(String params) throws IOException {
		return parse(new File(params));
	}

	public static MSFraggerParameter parse(File params) throws IOException {
		HashMap<String, MSFraggerParItem> valueMap = new HashMap<String, MSFraggerParItem>();

		int id = 1;
		String pre_anno = null;
		boolean blank = false;
		BufferedReader reader = new BufferedReader(new FileReader(params));
		String line = null;
		while ((line = reader.readLine()) != null) {
			if (line.startsWith("#")) {
				pre_anno = line;
				continue;
			}

			if (line.startsWith("variable_mod")) {
				continue;
			}

			if (line.trim().length() == 0) {
				blank = true;
			}

			int l1 = line.indexOf("=");
			if (l1 < 0) {
				continue;
			}

			String key = line.substring(0, l1).trim();
			String value = null;
			String annotation = null;
			int l2 = line.indexOf("#");
			if (l2 > l1) {
				value = line.substring(l1 + 1, l2).trim();
				annotation = line.substring(l2 + 1).trim();
			} else {
				value = line.substring(l1 + 1).trim();
			}

			MSFraggerParItem item = new MSFraggerParItem(id++, key, value, pre_anno, annotation, blank);
			valueMap.put(key, item);

			pre_anno = null;
			blank = false;
		}
		reader.close();

		MSFraggerParameter parameter = new MSFraggerParameter(valueMap);
		return parameter;
	}

	public static String getDefaultPath() {
		return params.getAbsolutePath();
	}

	public void setEnzyme(Enzyme enzyme, int termini, int miss) {
		valueMap.get("search_enzyme_name").setValue(enzyme.getName());
		valueMap.get("search_enzyme_cutafter").setValue(enzyme.getCleaveAt());
		valueMap.get("search_enzyme_butnotafter").setValue(enzyme.getNotCleaveAt());
		valueMap.get("num_enzyme_termini").setValue(String.valueOf(termini));
		valueMap.get("allowed_missed_cleavage").setValue(String.valueOf(miss));
	}

	public void setEnzyme(String enzyme, String cutAfter, String notCut, int termini, int miss) {
		valueMap.get("search_enzyme_name").setValue(enzyme);
		valueMap.get("search_enzyme_cutafter").setValue(cutAfter);
		valueMap.get("search_enzyme_butnotafter").setValue(notCut);
		valueMap.get("num_enzyme_termini").setValue(String.valueOf(termini));
		valueMap.get("allowed_missed_cleavage").setValue(String.valueOf(miss));
	}

	public void setDatabase(String database) {
		valueMap.get("database_name").setValue(database);
	}

	public void setThread(int threadCount) {
		valueMap.get("num_threads").setValue(String.valueOf(threadCount));
	}

	/**
	 * for the unit, 0=Daltons, 1=ppm
	 * 
	 * @param massTol
	 * @param unit
	 * @param trueTol
	 * @param trueUnit
	 */
	public void setPrecursorMassTol(double massTol, int unit, double trueTol, int trueUnit) {
		if (this.valueMap.containsKey("precursor_mass_tolerance")) {
			this.valueMap.get("precursor_mass_tolerance").setValue(String.valueOf(massTol));
		} else {
			if (this.valueMap.containsKey("precursor_mass_lower")) {
				this.valueMap.get("precursor_mass_lower").setValue(String.valueOf(-massTol));
			}
			if (this.valueMap.containsKey("precursor_mass_upper")) {
				this.valueMap.get("precursor_mass_upper").setValue(String.valueOf(massTol));
			}
		}

		this.valueMap.get("precursor_mass_units").setValue(String.valueOf(unit));
		this.valueMap.get("precursor_true_tolerance").setValue(String.valueOf(trueTol));
		this.valueMap.get("precursor_true_units").setValue(String.valueOf(trueUnit));
	}

	public void setFragmentMassTol(double massTol, int unit) {
		this.valueMap.get("fragment_mass_tolerance").setValue(String.valueOf(massTol));
		this.valueMap.get("fragment_mass_units").setValue(String.valueOf(unit));
	}

	public void export(String file) throws IOException {
		export(new File(file));
	}

	public void export(File file) throws IOException {
		PrintWriter writer = new PrintWriter(file);
		MSFraggerParItem[] items = this.valueMap.values().toArray(new MSFraggerParItem[this.valueMap.size()]);
		Arrays.sort(items);

		int id = 0;

		for (; id < items.length; id++) {
			writer.println(items[id]);
			if (items[id].key.equals("clip_nTerm_M")) {
				id++;
				break;
			}
		}

		String[] keys = variModMap.keySet().toArray(new String[variModMap.size()]);
		Arrays.sort(keys);
		for (int i = 0; i < keys.length; i++) {
			writer.println(keys[i] + " = " + variModMap.get(keys[i]));
		}

		for (; id < items.length; id++) {
			writer.println(items[id]);
		}

		writer.close();
	}

	private static class MSFraggerParItem implements Comparable<MSFraggerParItem> {

		private int id;
		private String key;
		private String value;
		private String pre_annotation;
		private String annotation;
		private boolean blank;

		MSFraggerParItem(int id, String key, String value, String pre_annotation, String annotation, boolean blank) {
			this.id = id;
			this.key = key;
			this.value = value;
			this.pre_annotation = pre_annotation;
			this.annotation = annotation;
			this.blank = blank;
		}

		void setValue(String value) {
			this.value = value;
		}

		public String toString() {
			StringBuilder sb = new StringBuilder();
			if (pre_annotation != null) {
				sb.append("\n").append(pre_annotation).append("\n");
			}
			if (blank) {
				sb.append("\n");
			}
			sb.append(key).append(" = ").append(value);
			if (annotation != null) {
				sb.append("\t").append(" # ").append(annotation);
			}
			return sb.toString();
		}

		@Override
		public int compareTo(MSFraggerParItem o) {
			// TODO Auto-generated method stub
			return this.id - o.id;
		}
	}

	public static void main(String[] args) throws IOException {
		// TODO Auto-generated method stub

		MSFraggerParameter par = MSFraggerParameter.parse("D:\\MetaLab\\Project\\Resources\\MSFragger\\open_fragger.params");
		System.out.println(par == null);
	}
}
