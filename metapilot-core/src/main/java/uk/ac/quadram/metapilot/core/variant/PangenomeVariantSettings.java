package uk.ac.quadram.metapilot.core.variant;

import java.util.Locale;

import org.json.JSONObject;

/**
 * Parameters of the pangenome variant check, shared by the two workflows that
 * can run it (MetaPilot open search and pFind open search). Both parameter
 * classes hold one instance and both ParameterIOs read and write it through
 * {@link #load(JSONObject)} / {@link #save(JSONObject)}, so the keys and their
 * defaults are defined once.
 * <p>
 * The tolerance for the MS2-undecidable test is not here: open search reuses
 * {@code substitution_mass_tolerance_da}, which already feeds
 * {@code SequenceVariantCaller}'s isobaric list, and pFind has its own
 * {@code undecidable_tolerance_da}.
 */
public final class PangenomeVariantSettings {

	public enum Mode {
		/** Runs when the UHGP-100 resource is present, silently skipped otherwise. */
		AUTO,
		/** Runs; a missing resource is reported as a warning, never as a failure. */
		ON,
		OFF;

		static Mode parse(String s) {
			if (s == null) {
				return AUTO;
			}
			switch (s.trim().toLowerCase(Locale.ROOT)) {
			case "on":
			case "true":
				return ON;
			case "off":
			case "false":
				return OFF;
			default:
				return AUTO;
			}
		}
	}

	public enum SpeciesScope {
		ALL, DETECTED, BOTH;

		static SpeciesScope parse(String s) {
			if (s == null) {
				return BOTH;
			}
			switch (s.trim().toLowerCase(Locale.ROOT)) {
			case "all":
				return ALL;
			case "detected":
				return DETECTED;
			default:
				return BOTH;
			}
		}

		public boolean reportsAll() {
			return this != DETECTED;
		}

		public boolean reportsDetected() {
			return this != ALL;
		}
	}

	public static final String KEY_MODE = "pangenome_variant_check";
	public static final String KEY_TRYPTIC_CONTEXT = "pangenome_require_tryptic_context";
	public static final String KEY_SPECIES_SCOPE = "pangenome_species_scope_report";
	public static final String KEY_CONTROL_SEED = "pangenome_control_seed";
	public static final String KEY_CONTROLS_PER_CANDIDATE = "pangenome_controls_per_candidate";
	/**
	 * The UHGP-100 resource: a prepared directory (holding {@code uhgp-100.faa.gz}),
	 * the {@code uhgp-100.faa.gz} itself, or the downloaded {@code uhgp-100.tar.gz},
	 * which is prepared once beside itself. Empty: look for
	 * {@code protein_catalogue/} under the selected MAG database version.
	 */
	public static final String KEY_RESOURCE_PATH = "uhgp100_path";
	/**
	 * MGnify's {@code uhgp-100.tsv} cluster membership (or {@code .tsv.gz}), or
	 * the directory holding it; it gives the genome counts. Empty: look beside
	 * {@code uhgp-100.faa.gz}, directly or in {@code uhgp-100/}.
	 */
	public static final String KEY_MEMBERSHIP_PATH = "uhgp100_membership_path";

	private Mode mode = Mode.AUTO;
	private boolean requireTrypticContext = true;
	private SpeciesScope speciesScope = SpeciesScope.BOTH;
	private long controlSeed = 1L;
	private int controlsPerCandidate = 1;
	private String resourcePath = "";
	private String membershipPath = "";

	public Mode getMode() {
		return mode;
	}

	public void setMode(Mode mode) {
		this.mode = mode == null ? Mode.AUTO : mode;
	}

	public boolean isRequireTrypticContext() {
		return requireTrypticContext;
	}

	public void setRequireTrypticContext(boolean requireTrypticContext) {
		this.requireTrypticContext = requireTrypticContext;
	}

	public SpeciesScope getSpeciesScope() {
		return speciesScope;
	}

	public void setSpeciesScope(SpeciesScope speciesScope) {
		this.speciesScope = speciesScope == null ? SpeciesScope.BOTH : speciesScope;
	}

	public long getControlSeed() {
		return controlSeed;
	}

	public void setControlSeed(long controlSeed) {
		this.controlSeed = controlSeed;
	}

	public int getControlsPerCandidate() {
		return controlsPerCandidate;
	}

	public void setControlsPerCandidate(int controlsPerCandidate) {
		this.controlsPerCandidate = Math.max(0, controlsPerCandidate);
	}

	public String getResourcePath() {
		return resourcePath;
	}

	public void setResourcePath(String resourcePath) {
		this.resourcePath = resourcePath == null ? "" : resourcePath.trim();
	}

	public String getMembershipPath() {
		return membershipPath;
	}

	public void setMembershipPath(String membershipPath) {
		this.membershipPath = membershipPath == null ? "" : membershipPath.trim();
	}

	/** Takes every value from {@code other}. */
	public void copyFrom(PangenomeVariantSettings other) {
		if (other == null) {
			return;
		}
		mode = other.mode;
		requireTrypticContext = other.requireTrypticContext;
		speciesScope = other.speciesScope;
		controlSeed = other.controlSeed;
		controlsPerCandidate = other.controlsPerCandidate;
		resourcePath = other.resourcePath;
		membershipPath = other.membershipPath;
	}

	/** Reads the keys from {@code additional_parameters}; absent keys keep their current value. */
	public void load(JSONObject extra) {
		if (extra == null) {
			return;
		}
		setMode(Mode.parse(extra.optString(KEY_MODE, mode.name())));
		setRequireTrypticContext(extra.optBoolean(KEY_TRYPTIC_CONTEXT, requireTrypticContext));
		setSpeciesScope(SpeciesScope.parse(extra.optString(KEY_SPECIES_SCOPE, speciesScope.name())));
		setControlSeed(extra.optLong(KEY_CONTROL_SEED, controlSeed));
		setControlsPerCandidate(extra.optInt(KEY_CONTROLS_PER_CANDIDATE, controlsPerCandidate));
		setResourcePath(extra.optString(KEY_RESOURCE_PATH, resourcePath));
		setMembershipPath(extra.optString(KEY_MEMBERSHIP_PATH, membershipPath));
	}

	/** Writes every key into {@code additional_parameters}. */
	public void save(JSONObject extra) {
		extra.put(KEY_MODE, mode.name().toLowerCase(Locale.ROOT));
		extra.put(KEY_TRYPTIC_CONTEXT, requireTrypticContext);
		extra.put(KEY_SPECIES_SCOPE, speciesScope.name().toLowerCase(Locale.ROOT));
		extra.put(KEY_CONTROL_SEED, controlSeed);
		extra.put(KEY_CONTROLS_PER_CANDIDATE, controlsPerCandidate);
		extra.put(KEY_RESOURCE_PATH, resourcePath);
		extra.put(KEY_MEMBERSHIP_PATH, membershipPath);
	}
}
