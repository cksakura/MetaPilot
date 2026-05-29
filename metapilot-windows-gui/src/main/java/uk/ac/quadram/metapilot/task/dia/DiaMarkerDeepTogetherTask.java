package uk.ac.quadram.metapilot.task.dia;

import javafx.scene.control.ProgressBar;

import uk.ac.quadram.metapilot.resources.MetaResourceParameters;
import uk.ac.quadram.metapilot.task.dia.par.MetaParameterDiaMarker;

/** Canonical Marker + Deep Together DIA workflow task. */
public class DiaMarkerDeepTogetherTask extends DiaMarkerDeepPerFileTask {

	public DiaMarkerDeepTogetherTask(MetaParameterDiaMarker markerParams, MetaResourceParameters mrp,
			ProgressBar bar1, ProgressBar bar2) {
		super(markerParams, mrp, bar1, bar2, DiaRefinementMode.DEEP_TOGETHER);
	}
}
