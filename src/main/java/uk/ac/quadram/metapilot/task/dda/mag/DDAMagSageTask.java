package uk.ac.quadram.metapilot.task.dda.mag;

import org.apache.logging.log4j.Logger;

import uk.ac.quadram.metapilot.gui.JavaFXProgressBar;
import uk.ac.quadram.metapilot.resources.MetaResourceParameters;
import uk.ac.quadram.metapilot.task.MetaPilotTask;
import uk.ac.quadram.metapilot.task.mag.MetaPilotMagTask;
import uk.ac.quadram.metapilot.task.par.MetaBaseParameter;

public class DDAMagSageTask extends MetaPilotMagTask {

	public DDAMagSageTask(MetaBaseParameter metaPar, MetaResourceParameters mrp, JavaFXProgressBar progressBar,
			JavaFXProgressBar secondBar) {
		super(metaPar, mrp, progressBar, secondBar);
	}

	@Override
	public void forceStop() {
		// TODO Auto-generated method stub
		
	}

	@Override
	protected String getTaskName() {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	protected Logger getLogger() {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	protected Boolean call() throws Exception {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	protected String getVersion() {
		// TODO Auto-generated method stub
		return null;
	}
	
	
}
