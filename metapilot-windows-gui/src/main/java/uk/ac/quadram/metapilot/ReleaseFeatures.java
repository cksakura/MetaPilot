package uk.ac.quadram.metapilot;

import uk.ac.quadram.metapilot.workflow.WorkflowType;

/**
 * Workflows that are built but not yet offered in the GUI.
 */
public final class ReleaseFeatures {

	/**
	 * The MetaPilot open search ({@link WorkflowType#DDA_MAG_OPENSEARCH}). Off hides its workflow
	 * card and result-cleaner entry and refuses its saved parameter files; starting the GUI with
	 * {@code -Dmetapilot.opensearch=true} turns it on.
	 */
	public static final boolean OPEN_SEARCH_WORKFLOW =
			!"false".equalsIgnoreCase(System.getProperty("metapilot.opensearch"));

	private ReleaseFeatures() {
	}

	public static boolean isAvailable(WorkflowType type) {
		return type != WorkflowType.DDA_MAG_OPENSEARCH || OPEN_SEARCH_WORKFLOW;
	}
}
