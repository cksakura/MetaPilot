package uk.ac.quadram.metapilot.task.dia;

import javafx.scene.control.ProgressBar;

import uk.ac.quadram.metapilot.resources.MetaResourceParameters;
import uk.ac.quadram.metapilot.task.dia.par.MetaParameterDia;

/**
 * DIA Normal workflow task (no ML model).
 * Extends MetaDiaTask for library-based identification and iterative search.
 */
public class MetaDiaNormalTask extends MetaDiaTask {

	private static final String TASK_NAME = "DIA Normal analysis";

	public MetaDiaNormalTask(MetaParameterDia metaPar, MetaResourceParameters mrp, ProgressBar bar1, ProgressBar bar2) {
		super(metaPar, mrp, bar1, bar2);
	}

	@Override
	protected String getTaskName() {
		return TASK_NAME;
	}
}
