package uk.ac.quadram.metapilot.task.dia;

import java.io.File;

import javafx.scene.control.ProgressBar;

import uk.ac.quadram.metapilot.core.runner.DiaMarkerDeepEntrapmentHeadlessRunner;
import uk.ac.quadram.metapilot.resources.MetaResourceParameters;
import uk.ac.quadram.metapilot.resources.MetaResourceParameters.ToolType;
import uk.ac.quadram.metapilot.task.dia.par.MetaParameterDiaMarker;
import uk.ac.quadram.metapilot.task.mag.MagDbItemResolver;

/** Canonical Marker + Deep Together DIA workflow task. */
public class DiaMarkerDeepTogetherTask extends DiaMarkerDeepPerFileTask {

	public DiaMarkerDeepTogetherTask(MetaParameterDiaMarker markerParams, MetaResourceParameters mrp,
			ProgressBar bar1, ProgressBar bar2) {
		super(markerParams, mrp, bar1, bar2, DiaRefinementMode.DEEP_TOGETHER);
	}

	@Override
	protected Boolean call() throws Exception {
		if (metaPar instanceof MetaParameterDiaMarker
				&& ((MetaParameterDiaMarker) metaPar).isEntrapmentEnabled()) {
			return callEntrapment();
		}
		return super.call();
	}

	private Boolean callEntrapment() throws Exception {
		MetaParameterDiaMarker markerParams = (MetaParameterDiaMarker) metaPar;
		logger.info(getTaskName() + ": starting Deep together paired-shuffle entrapment (seed "
				+ markerParams.getEntrapmentSeed() + ", ID digit " + markerParams.getEntrapmentIdDigit() + ")");
		updateMessage(getTaskName() + " (entrapment evaluation)");
		updateProgress(0, 100);
		updateSecondaryProgress(ProgressBar.INDETERMINATE_PROGRESS);

		if (!copyReportTemplatesInBackground()) {
			logger.error(getTaskName() + ": Failed to copy report templates");
			updateSecondaryProgress(0.0);
			return false;
		}

		updateMessage(getTaskName() + ": Converting spectra format...");
		updateProgress(10, 100);
		if (!spectraConvert()) {
			logger.error(getTaskName() + ": error in spectra format conversion");
			updateMessage(getTaskName() + ": Failed - Spectra conversion error");
			updateSecondaryProgress(0.0);
			return false;
		}
		updateProgress(25, 100);

		File diannExe = resolveGuiToolFile(ToolType.DIANN);
		if (diannExe == null || !diannExe.exists()) {
			logger.error(getTaskName() + ": DIA-NN executable not found");
			updateMessage(getTaskName() + ": Failed - DIA-NN executable not found");
			updateSecondaryProgress(0.0);
			return false;
		}

		try {
			File deepDir = DiaMarkerDeepEntrapmentHeadlessRunner.run(
					buildCoreRunnerParams(),
					diannExe,
					MagDbItemResolver.resolveConfigFile(null),
					createGuiProgressCallback(),
					inferCondaPrefixFromPython(resolveGuiToolFile(ToolType.PYTHON)),
					resolveGuiToolFile(ToolType.FUNCDB),
					null,
					markerParams.getEntrapmentSeed(),
					markerParams.getEntrapmentIdDigit());
			this.resultFolderFile = deepDir;
			this.metaPar.setDbSearchResultFile(deepDir);
			this.combinedFile = new File(deepDir, "combined");
			this.quan_pep_file = new File(this.combinedFile, "combined.pr_matrix.tsv");
			this.quan_pro_file = new File(this.combinedFile, "combined.pg_matrix.tsv");

			updateProgress(100, 100);
			updateMessage(getTaskName() + ": Entrapment evaluation completed ("
					+ DiaMarkerDeepEntrapmentHeadlessRunner.RESULT_FOLDER + ")");
			updateSecondaryProgress(1.0);
			return true;
		} catch (Exception e) {
			logger.error(getTaskName() + ": entrapment runner failed", e);
			updateMessage(getTaskName() + ": Failed - " + e.getMessage());
			updateSecondaryProgress(0.0);
			throw e;
		}
	}
}
