package uk.ac.quadram.metapilot.task.par;

import uk.ac.quadram.metapilot.task.dda.igc.par.DDAIgcMaxQuantParameterIO;
import uk.ac.quadram.metapilot.task.dda.igc.par.DDAIgcPFindParameterIO;
import uk.ac.quadram.metapilot.task.dda.mag.par.DDAMagFragPipeParameterIO;
import uk.ac.quadram.metapilot.task.dda.mag.par.DDAMagPfindParameterIO;
import uk.ac.quadram.metapilot.task.dda.mag.par.DDAMagSageParameterIO;
import uk.ac.quadram.metapilot.task.dia.par.DiaMarkerParameterIO;
import uk.ac.quadram.metapilot.task.dia.par.DiaParameterIO;
import uk.ac.quadram.metapilot.task.dia.par.DiaRTParameterIO;
import uk.ac.quadram.metapilot.task.dia.par.DiaTargetParameterIO;
import uk.ac.quadram.metapilot.workflow.WorkflowType;

/**
 * Registry for all ParameterIO classes. Call once before loading study JSON.
 */
public final class ParameterIORegistry {

    private ParameterIORegistry() {
    }

    public static void registerAllParameterIOs() {
        BaseParameterIO.registerParameterIO("DiaNNMAG", DiaParameterIO.class);
        BaseParameterIO.registerParameterIO("DiaNNMAGTG", DiaTargetParameterIO.class);
        BaseParameterIO.registerParameterIO("DiaNNMAGRT", DiaRTParameterIO.class);

        BaseParameterIO.registerParameterIO(WorkflowType.DIA_SEED.getAnalysisType(), DiaParameterIO.class);
        BaseParameterIO.registerParameterIO(WorkflowType.DIA_TARGET.getAnalysisType(), DiaTargetParameterIO.class);
        BaseParameterIO.registerParameterIO(WorkflowType.DIA_REALTIME.getAnalysisType(), DiaRTParameterIO.class);
        BaseParameterIO.registerParameterIO(WorkflowType.DIA_MARKER.getAnalysisType(), DiaMarkerParameterIO.class);

        BaseParameterIO.registerParameterIO(WorkflowType.DDA_MAG_PFIND.getAnalysisType(), DDAMagPfindParameterIO.class);
        BaseParameterIO.registerParameterIO(WorkflowType.DDA_MAG_FRAGPIPE.getAnalysisType(), DDAMagFragPipeParameterIO.class);
        BaseParameterIO.registerParameterIO(WorkflowType.DDA_MAG_SAGE.getAnalysisType(), DDAMagSageParameterIO.class);

        BaseParameterIO.registerParameterIO(WorkflowType.DDA_IGC_PFIND.getAnalysisType(), DDAIgcPFindParameterIO.class);
        BaseParameterIO.registerParameterIO(WorkflowType.DDA_IGC_MAXQUANT.getAnalysisType(), DDAIgcMaxQuantParameterIO.class);
        BaseParameterIO.registerParameterIO(WorkflowType.DDA_IGC_FRAGPIPE.getAnalysisType(), DDAMagFragPipeParameterIO.class);
        BaseParameterIO.registerParameterIO(WorkflowType.DDA_IGC_SAGE.getAnalysisType(), DDAMagSageParameterIO.class);
    }

    public static void printSupportedWorkflowTypes() {
        System.out.println("Supported workflow types:");
        for (String workflowType : BaseParameterIO.getSupportedWorkflowTypes()) {
            System.out.println("  - " + workflowType);
        }
    }
}
