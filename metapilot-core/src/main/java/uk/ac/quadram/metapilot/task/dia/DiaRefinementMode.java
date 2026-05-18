package uk.ac.quadram.metapilot.task.dia;

import java.util.Locale;

/**
 * Canonical DIA MAG refinement mode.
 * <p>
 * Persisted parameter files still use legacy values where needed:
 * FOCUS for {@link #CORE}, and DEEP for {@link #DEEP_TOGETHER}.
 */
public enum DiaRefinementMode {
	DIRECT("DIRECT", "Direct"),
	CORE("FOCUS", "Core / Focus"),
	COMPREHENSIVE("COMPREHENSIVE", "Comprehensive"),
	DEEP_TOGETHER("DEEP", "Deep together"),
	DEEP_PER_FILE("DEEP_PER_FILE", "Deep per file");

	private final String persistedValue;
	private final String displayLabel;

	DiaRefinementMode(String persistedValue, String displayLabel) {
		this.persistedValue = persistedValue;
		this.displayLabel = displayLabel;
	}

	public String getPersistedValue() {
		return persistedValue;
	}

	public String getDisplayLabel() {
		return displayLabel;
	}

	public boolean isDeep() {
		return this == DEEP_TOGETHER || this == DEEP_PER_FILE;
	}

	public boolean usesComprehensiveSelection() {
		return this == COMPREHENSIVE || this == DEEP_TOGETHER || this == DEEP_PER_FILE;
	}

	public static DiaRefinementMode fromPersistedValue(String value, DiaRefinementMode defaultMode) {
		if (value == null || value.isBlank()) {
			return defaultMode;
		}
		String normalized = value.trim().toUpperCase(Locale.ROOT).replace('-', '_').replace(' ', '_');
		switch (normalized) {
		case "DIRECT":
			return DIRECT;
		case "CORE":
		case "FOCUS":
			return CORE;
		case "COMPREHENSIVE":
			return COMPREHENSIVE;
		case "DEEP":
		case "DEEP_TOGETHER":
		case "DEEPTOGETHER":
			return DEEP_TOGETHER;
		case "DEEP_PER_FILE":
		case "DEEP_PER_RAW_FILE":
		case "DEEP_PERFILE":
			return DEEP_PER_FILE;
		default:
			return defaultMode;
		}
	}

	public static DiaRefinementMode fromMarkerMode(String markerMode) {
		return fromPersistedValue(markerMode, CORE);
	}

	public static DiaRefinementMode fromSeedMode(String seedMode) {
		return fromPersistedValue(seedMode, DIRECT);
	}

	public static DiaRefinementMode fromDisplayLabel(String label, DiaRefinementMode defaultMode) {
		if (label == null || label.isBlank()) {
			return defaultMode;
		}
		String normalized = label.trim().toLowerCase(Locale.ROOT);
		if (normalized.contains("direct")) {
			return DIRECT;
		}
		if (normalized.contains("per file") || normalized.contains("per raw file")) {
			return DEEP_PER_FILE;
		}
		if (normalized.contains("deep")) {
			return DEEP_TOGETHER;
		}
		if (normalized.contains("comprehensive")) {
			return COMPREHENSIVE;
		}
		if (normalized.contains("core") || normalized.contains("focus")) {
			return CORE;
		}
		return defaultMode;
	}
}
