package uk.ac.quadram.metapilot.core.runner;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import uk.ac.quadram.metapilot.core.dia.DiaDeepProteinMappingService;
import uk.ac.quadram.metapilot.core.workflow.ProgressCallback;
import uk.ac.quadram.metapilot.task.dia.DiaRefinementMode;
import uk.ac.quadram.metapilot.task.dia.par.MetaParameterDiaMarker;

import java.io.File;
import java.io.IOException;
import java.util.HashMap;
import java.util.HashSet;

/**
 * Local-side companion to {@link DiaMagReportHeadlessRunner}: reads the round parquet files
 * from a completed Deep Together / Deep Per File DIA Marker run (on the machine where they
 * already exist) and writes a single compact peptide-to-proteins TSV
 * ({@link DiaDeepProteinMappingService#writeMapTsv}).
 *
 * That one small file is all {@code diann-mag-report} needs on HPC to correctly remap deep
 * peptides to every protein they matched across rounds (see
 * {@link DiaMagReportHeadlessRunner} class docs for why {@code combined.parquet} alone is not
 * enough) - so only this file, not the whole {@code marker} round-parquet folder tree, needs
 * to be transferred.
 */
public final class DiaMagPeptideProteinMapExportRunner {

	private static final Logger LOGGER = LogManager.getLogger(DiaMagPeptideProteinMapExportRunner.class);
	private static final String TASK_NAME = "diann-mag-peptide-map";

	private DiaMagPeptideProteinMapExportRunner() {
	}

	public static final class Result {
		public final File outTsv;
		public final File roundRoot;
		public final int peptideCount;

		Result(File outTsv, File roundRoot, int peptideCount) {
			this.outTsv = outTsv;
			this.roundRoot = roundRoot;
			this.peptideCount = peptideCount;
		}
	}

	/**
	 * @param markerDeepRootOverride  optional override for the {@code marker_deep} folder
	 *                                (the one containing {@code marker/} and {@code combined/});
	 *                                when null, defaults to {@code <params.getResult()>/marker_deep} -
	 *                                the same location the local run itself wrote to, so this is
	 *                                normally run with no override, directly against the params
	 *                                file the local run used.
	 * @param outTsv                  where to write the compact peptide-protein map.
	 */
	public static Result run(MetaParameterDiaMarker params, File markerDeepRootOverride, File outTsv,
			ProgressCallback cb) throws Exception {
		if (cb == null)
			cb = new ProgressCallback.NoOp();
		if (params == null)
			throw new IllegalArgumentException("params is null");
		if (outTsv == null)
			throw new IllegalArgumentException("outTsv is null");

		DiaRefinementMode markerMode = params.getMarkerRefinementMode();
		if (markerMode != DiaRefinementMode.DEEP_TOGETHER && markerMode != DiaRefinementMode.DEEP_PER_FILE) {
			throw new IOException(TASK_NAME + ": only meaningful for Deep Together / Deep Per File marker modes "
					+ "(combined.parquet's own protein groups are already complete for Core/Focus/Comprehensive - "
					+ "no remap step exists or is needed for those modes). This params file has marker mode: "
					+ markerMode);
		}

		File markerDeepDir = markerDeepRootOverride != null ? markerDeepRootOverride
				: new File(params.getResult(), "marker_deep");
		String[] diaFiles = params.getMetadata() != null ? params.getMetadata().getRawFiles() : null;

		DiaDeepProteinMappingService.Request request;
		File checkDir;
		if (markerMode == DiaRefinementMode.DEEP_TOGETHER) {
			File markerDir = new File(markerDeepDir, "marker");
			checkDir = markerDir;
			request = DiaDeepProteinMappingService.Request.deepTogether(markerDir, 0, diaFiles, "KR");
		} else {
			checkDir = new File(markerDeepDir, "marker");
			request = DiaDeepProteinMappingService.Request.deepPerFile(markerDeepDir, diaFiles, "KR");
		}
		if (!checkDir.isDirectory()) {
			throw new IOException(TASK_NAME + ": round parquet folder not found: " + checkDir.getAbsolutePath()
					+ " (expected under <result>/marker_deep from the local run - pass --marker-root to override)");
		}

		cb.onMessage(TASK_NAME + ": reading round parquets under " + checkDir.getAbsolutePath() + " ...");
		HashMap<String, HashSet<String>> pepProMap = DiaDeepProteinMappingService.buildPeptideProteinMap(request,
				LOGGER, TASK_NAME);
		if (pepProMap.isEmpty()) {
			throw new IOException(TASK_NAME + ": no peptide-protein mappings found under "
					+ checkDir.getAbsolutePath());
		}

		cb.onMessage(TASK_NAME + ": writing " + pepProMap.size() + " peptide-protein mappings -> "
				+ outTsv.getAbsolutePath());
		DiaDeepProteinMappingService.writeMapTsv(outTsv, pepProMap);

		cb.onComplete(outTsv.getAbsolutePath());
		return new Result(outTsv, checkDir, pepProMap.size());
	}
}
