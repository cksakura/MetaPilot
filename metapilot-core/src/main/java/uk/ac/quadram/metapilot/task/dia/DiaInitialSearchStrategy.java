package uk.ac.quadram.metapilot.task.dia;

/**
 * Canonical DIA MAG initial evidence source.
 */
public enum DiaInitialSearchStrategy {
	MARKER("Marker"),
	SEED("Seed");

	private final String displayName;

	DiaInitialSearchStrategy(String displayName) {
		this.displayName = displayName;
	}

	public String getDisplayName() {
		return displayName;
	}
}
