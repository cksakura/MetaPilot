package uk.ac.quadram.metapilot.task.dia;

import javafx.scene.control.ProgressBar;

import uk.ac.quadram.metapilot.resources.MetaResourceParameters;
import uk.ac.quadram.metapilot.task.dia.par.MetaParameterDiaMarker;

/** Canonical Seed + Comprehensive DIA workflow task. */
public class DiaSeedComprehensiveTask extends DiaSeedCoreTask {

	public DiaSeedComprehensiveTask(MetaParameterDiaMarker markerParams, MetaResourceParameters mrp, ProgressBar bar1,
			ProgressBar bar2) {
		super(markerParams, mrp, bar1, bar2, DiaRefinementMode.COMPREHENSIVE);
	}
}
