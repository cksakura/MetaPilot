package uk.ac.quadram.metapilot.core.variant;

import org.json.JSONObject;

/**
 * Parameters of the open-search MS2 form ranking ({@link VariantFormRanking}),
 * stored under {@code additional_parameters} beside the pangenome check's own.
 * pFind has none: its substitution calls are trusted as-is.
 */
public final class VariantFormRankingSettings {

	public static final String KEY_ENABLED = "variant_ranking";
	/**
	 * Top two hypotheses within this many units of spectral-entropy similarity
	 * (0-1) are AMBIGUOUS. Provisional default; Kai sets the value.
	 */
	public static final String KEY_MARGIN = "variant_ranking_margin";
	/** Fragment matching tolerance; the bundled open-search parameters use 20 ppm. */
	public static final String KEY_FRAGMENT_TOLERANCE_PPM = "variant_ranking_fragment_tolerance_ppm";

	public static final double DEFAULT_MARGIN = 0.05;
	public static final double DEFAULT_FRAGMENT_TOLERANCE_PPM = 20.0;

	private boolean enabled = true;
	private double margin = DEFAULT_MARGIN;
	private double fragmentTolerancePpm = DEFAULT_FRAGMENT_TOLERANCE_PPM;

	public boolean isEnabled() {
		return enabled;
	}

	public void setEnabled(boolean enabled) {
		this.enabled = enabled;
	}

	public double getMargin() {
		return margin;
	}

	public void setMargin(double margin) {
		this.margin = Double.isFinite(margin) ? Math.max(0.0, margin) : DEFAULT_MARGIN;
	}

	public double getFragmentTolerancePpm() {
		return fragmentTolerancePpm;
	}

	public void setFragmentTolerancePpm(double ppm) {
		this.fragmentTolerancePpm = Double.isFinite(ppm) && ppm > 0 ? ppm : DEFAULT_FRAGMENT_TOLERANCE_PPM;
	}

	/** Reads the keys from {@code additional_parameters}; absent keys keep their current value. */
	public void load(JSONObject extra) {
		if (extra == null) {
			return;
		}
		setEnabled(extra.optBoolean(KEY_ENABLED, enabled));
		setMargin(extra.optDouble(KEY_MARGIN, margin));
		setFragmentTolerancePpm(extra.optDouble(KEY_FRAGMENT_TOLERANCE_PPM, fragmentTolerancePpm));
	}

	/** Writes every key into {@code additional_parameters}. */
	public void save(JSONObject extra) {
		extra.put(KEY_ENABLED, enabled);
		extra.put(KEY_MARGIN, margin);
		extra.put(KEY_FRAGMENT_TOLERANCE_PPM, fragmentTolerancePpm);
	}
}
