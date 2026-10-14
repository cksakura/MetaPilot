package uk.ac.quadram.metapilot.task.par;

import java.util.ArrayList;
import java.util.List;

import uk.ac.quadram.metapilot.core.mod.IsobaricTag;
import uk.ac.quadram.metapilot.dbSearch.MaxQuant.MaxquantModification;

/**
 * Resolves the isobaric-labelling settings a pFind workflow needs (tag,
 * per-channel reporter labels, isotope correction factors) from what the GUI
 * and the saved-parameter JSON actually provide.
 *
 * The pFind panels only hand over the quantification-mode text ("Isobaric"),
 * sometimes the selected tag name, the isobaric label table and the MetaData
 * label titles; they never call setIsobaricTag / setIsobaric / setIsoCorFactor.
 * Everything below is derived from those values so the search (fixed
 * modifications in pFind.cfg) and the reporter-ion quantification receive a
 * consistent tag.
 */
public final class IsobaricLabelSettings {

	/** Columns after the label name in an isobaric table row: -2, -1, +1, +2. */
	public static final int CORRECTION_COLUMNS = 4;

	private IsobaricLabelSettings() {
	}

	/** true when the quantification mode denotes isobaric labelling. */
	public static boolean isIsobaricMode(String quantitativeMode) {
		return quantitativeMode != null && quantitativeMode.toLowerCase().contains("isobaric");
	}

	/**
	 * Map the GUI radio-button texts ("Label-free", "Isobaric", "Isotopic") onto
	 * the {@link MetaConstants} values the tasks compare against. Unknown values
	 * are returned unchanged.
	 */
	public static String normalizeQuantitativeMode(String quantitativeMode) {
		if (quantitativeMode == null) {
			return null;
		}
		String m = quantitativeMode.trim().toLowerCase();
		if (m.contains("isobaric")) {
			return MetaConstants.isobaricLabel;
		}
		if (m.contains("isotopic")) {
			return MetaConstants.isotopicLabel;
		}
		if (m.contains("label-free") || m.contains("label free") || m.equals("lfq")) {
			return MetaConstants.labelFree;
		}
		return quantitativeMode;
	}

	/**
	 * Resolve the tag. Order: an explicitly set tag, the selected type name, the
	 * label names in the isobaric table, the MetaData label titles, the label
	 * types in the metadata table ("Isobaric - TMT10plex126C"). Returns null when
	 * the mode is not isobaric or nothing identifies a tag.
	 */
	public static IsobaricTag resolveTag(String quantitativeMode, IsobaricTag explicitTag, String isobaricType,
			List<String[]> isobaricTableData, MetaData metadata, List<String[]> metadataTableData) {

		if (!isIsobaricMode(quantitativeMode)) {
			return null;
		}
		if (explicitTag != null) {
			return explicitTag;
		}

		IsobaricTag tag = IsobaricTag.fromName(isobaricType);
		if (tag != null) {
			return tag;
		}

		tag = IsobaricTag.fromLabels(column(isobaricTableData, 0));
		if (tag != null) {
			return tag;
		}

		if (metadata != null) {
			tag = IsobaricTag.fromLabels(metadata.getLabelTitle());
			if (tag != null) {
				return tag;
			}
		}

		return IsobaricTag.fromLabels(column(metadataTableData, 1));
	}

	/**
	 * One reporter label per channel, in {@link IsobaricTag#getReportTitle()}
	 * order, carrying the neutral reporter mass as tmtMass.
	 */
	public static MaxquantModification[] buildReporterLabels(IsobaricTag tag) {
		if (tag == null) {
			return new MaxquantModification[0];
		}
		String[] titles = tag.getReportTitle();
		double[] reporter = tag.getReporterNeutralMass();
		String position = tag.isCysteineLabel() ? "anywhere" : "anyNterm";
		String[] sites = tag.isCysteineLabel() ? new String[] { "C" } : new String[] { "-" };

		MaxquantModification[] labels = new MaxquantModification[titles.length];
		for (int i = 0; i < titles.length; i++) {
			labels[i] = new MaxquantModification(titles[i], tag.getUnimodTitle() + " " + titles[i], position, "",
					sites, tag.getMass(), reporter[i], "");
		}
		return labels;
	}

	/**
	 * Isotope-impurity correction factors (percent), nPlex x 4 in
	 * {@link IsobaricTag#getReportTitle()} order. Rows are taken from the
	 * isobaric table ([label, -2, -1, +1, +2]) by label name; channels without a
	 * row, and blank or non-numeric cells, get 0 (no correction).
	 */
	public static double[][] buildCorrectionFactors(IsobaricTag tag, List<String[]> isobaricTableData) {
		if (tag == null) {
			return new double[0][0];
		}
		String[] titles = tag.getReportTitle();
		double[][] factors = new double[titles.length][CORRECTION_COLUMNS];
		if (isobaricTableData == null) {
			return factors;
		}
		for (String[] row : isobaricTableData) {
			if (row == null || row.length == 0 || row[0] == null) {
				continue;
			}
			String label = row[0].trim();
			for (int i = 0; i < titles.length; i++) {
				if (titles[i].equals(label)) {
					for (int j = 0; j < CORRECTION_COLUMNS; j++) {
						factors[i][j] = parse(row.length > j + 1 ? row[j + 1] : null);
					}
					break;
				}
			}
		}
		return factors;
	}

	/** true when the array has one row of {@link #CORRECTION_COLUMNS} per channel of the tag. */
	public static boolean correctionFactorsFit(IsobaricTag tag, double[][] factors) {
		if (tag == null || factors == null || factors.length != tag.getReportTitle().length) {
			return false;
		}
		for (double[] row : factors) {
			if (row == null || row.length != CORRECTION_COLUMNS) {
				return false;
			}
		}
		return true;
	}

	/** true when the labels are one per channel of the tag, in report-title order. */
	public static boolean reporterLabelsFit(IsobaricTag tag, MaxquantModification[] labels) {
		if (tag == null || labels == null) {
			return false;
		}
		String[] titles = tag.getReportTitle();
		if (labels.length != titles.length) {
			return false;
		}
		for (int i = 0; i < titles.length; i++) {
			if (labels[i] == null || labels[i].getTmtMass() <= 0) {
				return false;
			}
		}
		return true;
	}

	private static String[] column(List<String[]> rows, int col) {
		if (rows == null || rows.isEmpty()) {
			return null;
		}
		List<String> values = new ArrayList<String>();
		for (String[] row : rows) {
			if (row != null && row.length > col && row[col] != null && !row[col].trim().isEmpty()) {
				values.add(row[col]);
			}
		}
		return values.isEmpty() ? null : values.toArray(new String[0]);
	}

	private static double parse(String s) {
		if (s == null || s.trim().isEmpty()) {
			return 0.0;
		}
		try {
			return Double.parseDouble(s.trim());
		} catch (NumberFormatException e) {
			return 0.0;
		}
	}
}
