package uk.ac.quadram.metapilot.task.par;

import uk.ac.quadram.metapilot.MetaPilotApp;
import uk.ac.quadram.metapilot.task.MetaLabWorkflowType;
import uk.ac.quadram.metapilot.task.dia.par.DiaParameterIO;
import uk.ac.quadram.metapilot.task.dia.par.DiaTargetParameterIO;
import uk.ac.quadram.metapilot.task.dia.par.DiaRTParameterIO;
import uk.ac.quadram.metapilot.task.dda.igc.par.DDAIgcMaxQuantParameterIO;
import uk.ac.quadram.metapilot.task.dda.igc.par.DDAIgcPFindParameterIO;
import uk.ac.quadram.metapilot.task.dda.mag.par.DDAMagPfindParameterIO;
import uk.ac.quadram.metapilot.task.dda.mag.par.DDAMagFragPipeParameterIO;
import uk.ac.quadram.metapilot.task.dda.mag.par.DDAMagSageParameterIO;

/**
 * Registry for all ParameterIO classes This class should be called during
 * application initialization to register all available workflow types with
 * their corresponding ParameterIO classes.
 */
public class ParameterIORegistry {

	/**
	 * Register all available ParameterIO classes This should be called once during
	 * application startup
	 */
	public static void registerAllParameterIOs() {
		// Register DIA workflow types using MetaLabWorkflowType enum
		BaseParameterIO.registerParameterIO(MetaPilotApp.WorkflowType.DIA_NORMAL.getAnalysisType(),
				DiaParameterIO.class);
		BaseParameterIO.registerParameterIO(MetaPilotApp.WorkflowType.DIA_TARGET.getAnalysisType(),
				DiaTargetParameterIO.class);
		BaseParameterIO.registerParameterIO(MetaPilotApp.WorkflowType.DIA_REALTIME.getAnalysisType(),
				DiaRTParameterIO.class);

		// Register DDA MAG workflow types
		BaseParameterIO.registerParameterIO(MetaPilotApp.WorkflowType.DDA_MAG_PFIND.getAnalysisType(),
				DDAMagPfindParameterIO.class);
		BaseParameterIO.registerParameterIO(MetaPilotApp.WorkflowType.DDA_MAG_FRAGPIPE.getAnalysisType(),
				DDAMagFragPipeParameterIO.class);
		BaseParameterIO.registerParameterIO(MetaPilotApp.WorkflowType.DDA_MAG_SAGE.getAnalysisType(),
				DDAMagSageParameterIO.class);

		// Register DDA IGC workflow types
		BaseParameterIO.registerParameterIO(MetaPilotApp.WorkflowType.DDA_IGC_PFIND.getAnalysisType(),
				DDAIgcPFindParameterIO.class);
		BaseParameterIO.registerParameterIO(MetaPilotApp.WorkflowType.DDA_IGC_MAXQUANT.getAnalysisType(),
				DDAIgcMaxQuantParameterIO.class);
		BaseParameterIO.registerParameterIO(MetaPilotApp.WorkflowType.DDA_IGC_FRAGPIPE.getAnalysisType(),
				DDAMagFragPipeParameterIO.class);
		BaseParameterIO.registerParameterIO(MetaPilotApp.WorkflowType.DDA_IGC_SAGE.getAnalysisType(),
				DDAMagSageParameterIO.class);
	}

	/**
	 * Get a list of all registered workflow types for CLI help
	 */
	public static void printSupportedWorkflowTypes() {
		System.out.println("Supported workflow types:");
		for (String workflowType : BaseParameterIO.getSupportedWorkflowTypes()) {
			System.out.println("  - " + workflowType);
		}
	}

	/**
	 * Get the workflow type enum from a string name
	 */
	public static MetaLabWorkflowType getWorkflowTypeFromName(String workflowName) {
		try {
			return MetaLabWorkflowType.valueOf(workflowName);
		} catch (IllegalArgumentException e) {
			return null;
		}
	}
}
