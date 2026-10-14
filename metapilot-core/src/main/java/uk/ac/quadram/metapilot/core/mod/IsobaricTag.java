/**
 * 
 */
package uk.ac.quadram.metapilot.core.mod;

import java.util.Arrays;
import java.util.LinkedHashMap;

import uk.ac.quadram.metapilot.core.aminoacid.AminoAcidProperty;

/**
 * @author Kai Cheng
 *
 */
public enum IsobaricTag {

	itraq4("iTRAQ4plex", "iTRAQ (Representative mass and accurate mass for 116 & 117)", 144.102063, 4),

	itraq8("iTRAQ8plex", "iTRAQ8plex (Representative mass and accurate mass for 113, 114, 116 & 117)", 304.205360, 8),

	tmt2("TMT2plex", "Duplex Tandem Mass Tag", 225.155833, 2),

	tmt6("TMT6plex", "Sixplex Tandem Mass Tag", 229.162932, 6),

	tmt8("TMT8plex", "TMT8plex", 229.162932, 8),

	tmt10("TMT10plex", "TMT10plex", 229.162932, 10),

	tmt11("TMT11plex", "TMT11plex", 229.162932, 11),

	iodo6("iodoTMT6plex", "Sixplex iodoacetyl Tandem Mass Tag", 329.226595, 6),

	tmtpro16("TMTpro16plex", "TMTpro 16plex", 304.207146, 16),

	tmtpro18("TMTpro18plex", "TMTpro 18plex", 304.207146, 18);

	private String name;
	private String unimodTitle;
	private double mass;
	private int nPlex;

	IsobaricTag(String name, String unimodTitle, double mass, int nPlex) {
		this.name = name;
		this.unimodTitle = unimodTitle;
		this.mass = mass;
		this.nPlex = nPlex;
	}

	public String getName() {
		return name;
	}

	public String getUnimodTitle() {
		return unimodTitle;
	}

	public double getMass() {
		return mass;
	}

	public int getnPlex() {
		return nPlex;
	}

	public String[] getReportTitle() {
		switch (this) {
		case itraq4: {
			String[] ions = new String[4];
			String[] itaq = itraqIons();
			System.arraycopy(itaq, 1, ions, 0, ions.length);
			for (int i = 0; i < ions.length; i++) {
				ions[i] = this.getName() + ions[i];
			}
			return ions;
		}
		case itraq8: {
			String[] ions = itraqIons();
			for (int i = 0; i < ions.length; i++) {
				ions[i] = this.getName() + ions[i];
			}
			return ions;
		}
		case tmt2: {
			String[] ions = new String[2];
			String[] tmt = tmtIons();
			System.arraycopy(tmt, 1, ions, 0, ions.length);
			for (int i = 0; i < ions.length; i++) {
				ions[i] = this.getName() + ions[i];
			}
			return ions;
		}
		case tmt6: {
			String[] ions = tmtIons();
			for (int i = 0; i < ions.length; i++) {
				ions[i] = this.getName() + ions[i];
			}
			return ions;
		}
		case tmt8: {
			String[] tmt = tmtIonsNC();
			String[] ions = new String[] { tmt[0], tmt[1], tmt[2], tmt[4], tmt[5], tmt[6], tmt[8], tmt[9] };
			for (int i = 0; i < ions.length; i++) {
				ions[i] = this.getName() + ions[i];
			}
			return ions;
		}
		case tmt10: {
			String[] ions = new String[10];
			String[] tmt = tmtIonsNC();
			System.arraycopy(tmt, 0, ions, 0, ions.length);
			for (int i = 0; i < ions.length; i++) {
				ions[i] = this.getName() + ions[i];
			}
			return ions;
		}
		case tmt11: {
			String[] ions = tmtIonsNC();
			for (int i = 0; i < ions.length - 1; i++) {
				ions[i] = tmt10.getName() + ions[i];
			}
			ions[ions.length - 1] = this.getName() + ions[ions.length - 1];
			return ions;
		}
		case iodo6: {
			String[] ions = iodoTmtIons();
			for (int i = 0; i < ions.length; i++) {
				ions[i] = this.getName() + ions[i];
			}
			return ions;
		}
		case tmtpro16: {
			String[] ions = new String[16];
			String[] tmt = tmtProIonsNC();
			System.arraycopy(tmt, 0, ions, 0, ions.length);
			for (int i = 0; i < ions.length; i++) {
				ions[i] = this.getName() + ions[i];
			}
			return ions;
		}
		case tmtpro18: {
			String[] ions = tmtProIonsNC();
			for (int i = 0; i < ions.length; i++) {
				ions[i] = this.getName() + ions[i];
			}
			return ions;
		}
		default: {
			return new String[] {};
		}
		}
	}

	public static IsobaricTag getTagFromReportTitle(String[] reportTitle) {
		L: for (IsobaricTag tag : IsobaricTag.values()) {
			String[] tagTitle = tag.getReportTitle();
			for (int i = 0; i < reportTitle.length; i++) {
				boolean find = false;
				for (int j = 0; j < tagTitle.length; j++) {
					if (reportTitle[i].equals(tagTitle[j])) {
						find = true;
						break;
					}
				}
				if (!find) {
					continue L;
				}
			}
			return tag;
		}
		return null;
	}

	private static String[] itraqIons() {
		String[] ions = new String[8];
		for (int i = 0; i < 7; i++) {
			ions[i] = String.valueOf(i + 113);
		}
		ions[7] = "121";
		return ions;
	}

	private static String[] tmtIons() {
		String[] ions = new String[6];
		for (int i = 0; i < 6; i++) {
			ions[i] = String.valueOf(i + 126);
		}
		return ions;
	}

	private static String[] tmtIonsNC() {
		String[] ions = new String[11];
		ions[0] = "126C";
		for (int i = 0; i < 5; i++) {
			ions[i * 2 + 1] = (i + 127) + "N";
			ions[i * 2 + 2] = (i + 127) + "C";
		}
		return ions;
	}

	private static String[] tmtProIonsNC() {
		String[] ions = new String[18];
		ions[0] = "126C";
		for (int i = 0; i < 8; i++) {
			ions[i * 2 + 1] = (i + 127) + "N";
			ions[i * 2 + 2] = (i + 127) + "C";
		}
		ions[ions.length - 1] = "135N";
		return ions;
	}

	private static String[] iodoTmtIons() {
		String[] ions = new String[6];
		for (int i = 0; i < 6; i++) {
			ions[i] = "Cys" + (i + 126);
		}
		return ions;
	}

	public static IsobaricTag[] getMaxQuantTags() {
		return new IsobaricTag[] { itraq4, itraq8, tmt2, tmt6, tmt8, tmt10, tmt11, iodo6 };
	}

	public static IsobaricTag[] getPFindTags() {
		return new IsobaricTag[] { itraq4, itraq8, tmt2, tmt6, tmt8, tmt10, tmt11, iodo6, tmtpro16, tmtpro18 };
	}
	
	public static IsobaricTag[] getFragpipeTags() {
		return new IsobaricTag[] { itraq4, itraq8, tmt6, tmt10, tmt11, tmtpro16, tmtpro18 };
	}

	// ------------------------------------------------------------------
	// Reporter ions
	// ------------------------------------------------------------------

	/**
	 * Reporter ion m/z (singly charged) of the TMT/TMTpro channels, in the order
	 * of {@link #tmtProIonsNC()}: 126C, 127N, 127C, 128N, 128C, 129N, 129C, 130N,
	 * 130C, 131N, 131C, 132N, 132C, 133N, 133C, 134N, 134C, 135N. The first 11
	 * entries are also the order of {@link #tmtIonsNC()}.
	 */
	private static final double[] TMT_NC_REPORTER_MZ = { 126.127726, 127.124761, 127.131081, 128.128116,
			128.134436, 129.131471, 129.137790, 130.134825, 130.141145, 131.138180, 131.144500, 132.141535,
			132.147855, 133.144890, 133.151210, 134.148245, 134.154565, 135.151600 };

	/** Reporter ion m/z of the iTRAQ channels 113, 114, 115, 116, 117, 118, 119, 121. */
	private static final double[] ITRAQ_REPORTER_MZ = { 113.107873, 114.111228, 115.108263, 116.111618, 117.114973,
			118.112008, 119.115363, 121.122072 };

	/**
	 * Reporter ion m/z (singly charged) of every channel, in the same order as
	 * {@link #getReportTitle()}.
	 */
	public double[] getReporterIonMz() {
		switch (this) {
		case itraq4:
			return Arrays.copyOfRange(ITRAQ_REPORTER_MZ, 1, 5);
		case itraq8:
			return ITRAQ_REPORTER_MZ.clone();
		case tmt2:
			// TMT duplex = 126 and 127C
			return new double[] { TMT_NC_REPORTER_MZ[0], TMT_NC_REPORTER_MZ[2] };
		case tmt6:
		case iodo6:
			// sixplex = 126, 127N, 128C, 129N, 130C, 131N
			return new double[] { TMT_NC_REPORTER_MZ[0], TMT_NC_REPORTER_MZ[1], TMT_NC_REPORTER_MZ[4],
					TMT_NC_REPORTER_MZ[5], TMT_NC_REPORTER_MZ[8], TMT_NC_REPORTER_MZ[9] };
		case tmt8:
			// same channel selection as getReportTitle(): 126C 127N 127C 128C 129N 129C 130C 131N
			return new double[] { TMT_NC_REPORTER_MZ[0], TMT_NC_REPORTER_MZ[1], TMT_NC_REPORTER_MZ[2],
					TMT_NC_REPORTER_MZ[4], TMT_NC_REPORTER_MZ[5], TMT_NC_REPORTER_MZ[6], TMT_NC_REPORTER_MZ[8],
					TMT_NC_REPORTER_MZ[9] };
		case tmt10:
			return Arrays.copyOf(TMT_NC_REPORTER_MZ, 10);
		case tmt11:
			return Arrays.copyOf(TMT_NC_REPORTER_MZ, 11);
		case tmtpro16:
			return Arrays.copyOf(TMT_NC_REPORTER_MZ, 16);
		case tmtpro18:
			return Arrays.copyOf(TMT_NC_REPORTER_MZ, 18);
		default:
			return new double[0];
		}
	}

	/**
	 * Neutral reporter masses (reporter m/z minus a proton), in the same order as
	 * {@link #getReportTitle()}. This is the value carried as "tmtMass" on the
	 * isobaric label objects; the reporter extraction adds the proton back.
	 */
	public double[] getReporterNeutralMass() {
		double[] mz = getReporterIonMz();
		double[] mass = new double[mz.length];
		for (int i = 0; i < mz.length; i++) {
			mass[i] = mz[i] - AminoAcidProperty.PROTON_W;
		}
		return mass;
	}

	// ------------------------------------------------------------------
	// pFind
	// ------------------------------------------------------------------

	/** true when the tag labels cysteine instead of the peptide amines. */
	public boolean isCysteineLabel() {
		return this == iodo6;
	}

	/**
	 * pFind modification names (modification.ini keys) written as fixed
	 * modifications for this tag.
	 */
	public String[] getPFindFixedMods() {
		switch (this) {
		case itraq4:
			return new String[] { "iTRAQ4plex[AnyN-term]", "iTRAQ4plex[K]" };
		case itraq8:
			return new String[] { "iTRAQ8plex[AnyN-term]", "iTRAQ8plex[K]" };
		case tmt2:
			return new String[] { "TMT2plex[AnyN-term]", "TMT2plex[K]" };
		case tmt6:
		case tmt8:
		case tmt10:
		case tmt11:
			return new String[] { "TMT6plex[AnyN-term]", "TMT6plex[K]" };
		case iodo6:
			return new String[] { "iodoTMT6plex[C]" };
		case tmtpro16:
		case tmtpro18:
			return new String[] { "TMTpro[AnyN-term]", "TMTpro[K]" };
		default:
			return new String[0];
		}
	}

	/**
	 * modification.ini definition values ("sites position mono avg 0
	 * composition") keyed by every name returned by {@link #getPFindFixedMods()}.
	 * Used to add missing entries to a copy of pFind's modification.ini: the
	 * pFind 3 ini has no TMTpro entries and no iTRAQ4plex[AnyN-term]/[K] entries.
	 */
	public LinkedHashMap<String, String> getPFindModDefinitions() {
		LinkedHashMap<String, String> map = new LinkedHashMap<String, String>();
		String nterm = "ABCDEFGHIJKLMNOPQRSTUVWXYZ PEP_N ";
		switch (this) {
		case itraq4: {
			String def = "144.102063 144.1544 0 H(12)C(4)13C(3)N(1)15N(1)O(1)";
			map.put("iTRAQ4plex[AnyN-term]", nterm + def);
			map.put("iTRAQ4plex[K]", "K NORMAL " + def);
			break;
		}
		case itraq8: {
			String def = "304.20536 304.3074 0 H(24)C(7)13C(7)N(3)15N(1)O(3)";
			map.put("iTRAQ8plex[AnyN-term]", nterm + def);
			map.put("iTRAQ8plex[K]", "K NORMAL " + def);
			break;
		}
		case tmt2: {
			String def = "225.155833 225.2921 0 H(20)C(11)13C(1)N(2)O(2)";
			map.put("TMT2plex[AnyN-term]", nterm + def);
			map.put("TMT2plex[K]", "K NORMAL " + def);
			break;
		}
		case tmt6:
		case tmt8:
		case tmt10:
		case tmt11: {
			String def = "229.162932 229.2634 0 H(20)C(8)13C(4)N(1)15N(1)O(2)";
			map.put("TMT6plex[AnyN-term]", nterm + def);
			map.put("TMT6plex[K]", "K NORMAL " + def);
			break;
		}
		case iodo6:
			map.put("iodoTMT6plex[C]", "C NORMAL 329.226595 329.3825 0 H(28)C(12)13C(4)N(3)15N(1)O(3)");
			break;
		case tmtpro16:
		case tmtpro18: {
			// Unimod #2016 TMTpro: H(25) C(8) 13C(7) N 15N(2) O(3)
			String def = "304.207146 304.3127 0 H(25)C(8)13C(7)N(1)15N(2)O(3)";
			map.put("TMTpro[AnyN-term]", nterm + def);
			map.put("TMTpro[K]", "K NORMAL " + def);
			break;
		}
		default:
			break;
		}
		return map;
	}

	// ------------------------------------------------------------------
	// Resolution from GUI / saved-parameter values
	// ------------------------------------------------------------------

	/**
	 * Find a tag by {@link #getName()} or enum constant name, ignoring case.
	 * Returns null for anything else, e.g. the quantification-mode text
	 * "Isobaric".
	 */
	public static IsobaricTag fromName(String name) {
		if (name == null) {
			return null;
		}
		String n = name.trim();
		if (n.isEmpty()) {
			return null;
		}
		for (IsobaricTag tag : IsobaricTag.values()) {
			if (tag.getName().equalsIgnoreCase(n) || tag.name().equalsIgnoreCase(n)) {
				return tag;
			}
		}
		return null;
	}

	/**
	 * Find the tag whose channel labels contain all of the given labels (e.g.
	 * "TMT10plex126C", "TMT10plex127N", ...). Null/blank entries and the
	 * "Isobaric - " prefix used in the metadata table are ignored. Returns null
	 * when no label is given or no tag matches.
	 */
	public static IsobaricTag fromLabels(String[] labels) {
		if (labels == null) {
			return null;
		}
		java.util.ArrayList<String> list = new java.util.ArrayList<String>();
		for (String label : labels) {
			if (label == null) {
				continue;
			}
			String l = label.trim();
			if (l.startsWith("Isobaric - ")) {
				l = l.substring("Isobaric - ".length()).trim();
			}
			if (!l.isEmpty()) {
				list.add(l);
			}
		}
		if (list.isEmpty()) {
			return null;
		}
		return getTagFromReportTitle(list.toArray(new String[0]));
	}
}
