package uk.ac.quadram.metapilot.task.dia.par;

import java.io.File;

import uk.ac.quadram.metapilot.task.mag.MagDbItem;
import uk.ac.quadram.metapilot.task.par.MetaBaseParameter;

public class MetaParameterDiaRT extends MetaParameterDia {

	private boolean combineResult;
	private boolean stopNow;
	private boolean stopComplete;
	private int taskCount;
	private File monitorFolder;
	private File resultFolder;
	/**
	 * Selected real-time workflow mode.
	 * <p>
	 * Expected values (string codes, interpreted by UI and monitor task):
	 * <ul>
	 *   <li>"MARKER_FOCUS"</li>
	 *   <li>"MARKER_COMPREHENSIVE"</li>
	 *   <li>"MARKER_DEEP"</li>
	 *   <li>"SEED_DIRECT"</li>
	 *   <li>"SEED_FOCUS"</li>
	 *   <li>"SEED_COMPREHENSIVE"</li>
	 *   <li>"SEED_DEEP"</li>
	 * </ul>
	 * Implementations should treat unknown or null as a sensible default (e.g. SEED_DIRECT).
	 */
	private String rtWorkflowMode;

	public MetaParameterDiaRT(MetaBaseParameter metaPar, MagDbItem magDbItem, boolean librarySearch,
			boolean selfModelSearch, String[] library, boolean useModel, String[] models) {
		super(metaPar, magDbItem, librarySearch, selfModelSearch, library, useModel, models);
		// TODO Auto-generated constructor stub
	}

	public boolean isCombineResult() {
		return combineResult;
	}

	public void setCombineResult(boolean combineResult) {
		this.combineResult = combineResult;
	}

	public boolean isStopNow() {
		return stopNow;
	}

	public void setStopNow(boolean stopNow) {
		this.stopNow = stopNow;
	}

	public boolean isStopComplete() {
		return stopComplete;
	}

	public void setStopComplete(boolean stopComplete) {
		this.stopComplete = stopComplete;
	}

	public int getTaskCount() {
		return taskCount;
	}

	public void setTaskCount(int taskCount) {
		this.taskCount = taskCount;
	}

	public File getMonitorFolder() {
		return monitorFolder;
	}

	public void setMonitorFolder(File monitorFolder) {
		this.monitorFolder = monitorFolder;
	}

	public File getResultFolder() {
		return resultFolder;
	}

	public void setResultFolder(File resultFolder) {
		this.resultFolder = resultFolder;
	}

	/**
	 * Returns the selected real-time workflow mode code.
	 */
	public String getRtWorkflowMode() {
		return rtWorkflowMode;
	}

	/**
	 * Sets the selected real-time workflow mode code.
	 *
	 * @param rtWorkflowMode one of the documented mode strings, or null to use a default.
	 */
	public void setRtWorkflowMode(String rtWorkflowMode) {
		this.rtWorkflowMode = rtWorkflowMode;
	}

}
