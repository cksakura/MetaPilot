package uk.ac.quadram.metapilot.task.io.qpx;

import java.io.File;
import java.util.Arrays;
import java.util.Collection;
import java.util.function.Consumer;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import uk.ac.quadram.metapilot.core.workflow.ProgressCallback;
import uk.ac.quadram.metapilot.task.io.pep.MetaPSM;

/**
 * The one call every DDA workflow makes to publish its PSM layer as QPX: writes
 * {@value QpxPsmParquetWriter#DEFAULT_FILE_NAME} beside {@code final_peptides.tsv}.
 * <p>
 * Never throws: the PSM parquet is an additional deliverable, so a failure is logged
 * and reported through the callback, and the workflow's own results stand.
 */
public final class QpxPsmExport {

	private static final Logger LOGGER = LogManager.getLogger(QpxPsmExport.class);

	private QpxPsmExport() {
	}

	public static File writeQuietly(File resultDir, Collection<? extends MetaPSM> psms, ProgressCallback cb) {
		return writeQuietly(resultDir, psms, cb == null ? null : (Consumer<String>) cb::onMessage);
	}

	public static File writeQuietly(File resultDir, MetaPSM[] psms, ProgressCallback cb) {
		return writeQuietly(resultDir, psms == null ? null : Arrays.asList(psms), cb);
	}

	/** Variant for GUI tasks that log through their own sink rather than a ProgressCallback. */
	public static File writeQuietly(File resultDir, Collection<? extends MetaPSM> psms, Consumer<String> log) {
		if (resultDir == null || psms == null) {
			return null;
		}
		File out = new File(resultDir, QpxPsmParquetWriter.DEFAULT_FILE_NAME);
		try {
			QpxPsmParquetWriter.Summary summary = QpxPsmParquetWriter.write(out, psms);
			String msg = "QPX PSM export: " + summary + " -> " + out.getAbsolutePath();
			LOGGER.info(msg);
			if (log != null) {
				log.accept(msg);
			}
			return out;
		} catch (Throwable t) {
			String msg = "QPX PSM export skipped (" + out.getAbsolutePath() + "): " + t;
			LOGGER.warn(msg, t);
			if (log != null) {
				log.accept(msg);
			}
			return null;
		}
	}
}
