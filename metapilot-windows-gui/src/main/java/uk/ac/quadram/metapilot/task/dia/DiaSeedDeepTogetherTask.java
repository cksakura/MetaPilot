package uk.ac.quadram.metapilot.task.dia;

import javafx.scene.control.ProgressBar;

import uk.ac.quadram.metapilot.resources.MetaResourceParameters;
import uk.ac.quadram.metapilot.task.dia.par.MetaParameterDiaMarker;

/** Canonical Seed + Deep Together DIA workflow task. */
public class DiaSeedDeepTogetherTask extends DiaSeedDeepPerFileTask {

	public DiaSeedDeepTogetherTask(MetaParameterDiaMarker markerParams, MetaResourceParameters mrp, ProgressBar bar1,
			ProgressBar bar2) {
		super(markerParams, mrp, bar1, bar2, DiaRefinementMode.DEEP_TOGETHER);
	}
}
