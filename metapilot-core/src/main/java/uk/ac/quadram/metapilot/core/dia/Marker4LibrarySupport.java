package uk.ac.quadram.metapilot.core.dia;

import java.io.File;
import java.io.PrintWriter;
import java.util.HashSet;

import org.apache.logging.log4j.Logger;

import uk.ac.quadram.metapilot.core.diann.DiaNNHeadlessRunner;
import uk.ac.quadram.metapilot.core.workflow.ProgressCallback;
import uk.ac.quadram.metapilot.dbSearch.diann.DiaNNParquetReader;
import uk.ac.quadram.metapilot.dbSearch.diann.DiannParameter;
import uk.ac.quadram.metapilot.task.dia.par.MetaParameterDia;
import uk.ac.quadram.metapilot.task.mag.MagDbItem;

/**
 * Shared helper for maintaining Marker4 library metadata under {@code libraries/Marker4/}.
 */
public final class Marker4LibrarySupport {

	private static final String ENZYME = "Trypsin/P";
	private static final int MISS = 1;
	private static final boolean MATCH_TO_MAG = true;

	@FunctionalInterface
	public interface SpeclibToParquetConverter {
		boolean convert(DiannParameter diannPar, File speclib, File parquet, File libFolder,
				ProgressCallback callback, Logger logger) throws Exception;
	}

	private Marker4LibrarySupport() {
	}

	private static DiannParameter buildDiannParForLibMetadata(MetaParameterDia diaPar, boolean speclibOnly) {
		DiannParameter diannPar = new DiannParameter();
		diannPar.setThreads(diaPar.getThreadCount());
		diannPar.setEnzyme(ENZYME);
		diannPar.setMissed_cleavages(MISS);
		diannPar.setM_excision(diaPar.isMExcision());
		diannPar.setCarbami_C(diaPar.isCarbamiC());
		diannPar.setOxidation_M(diaPar.isOxidationM());
		diannPar.setAcety_proN(diaPar.isAcetyProN());
		diannPar.setPhospho_STY(diaPar.isPhosphoSTY());
		diannPar.setK_GG(diaPar.isKGG());
		if (speclibOnly) {
			diannPar.setMatrices(false);
			diannPar.setReanalyse(true);
		}
		return diannPar;
	}

	/**
	 * Ensures {@code Marker4.parquet}, {@code Marker4_info.txt}, and in-memory registration exist for the canonical
	 * speclib (parquet + info + {@link MagDbItem#addPepLib}).
	 */
	public static boolean ensureMarker4LibraryMetadata(MagDbItem magDb, MetaParameterDia diaPar,
			ProgressCallback callback, Logger logger, boolean speclibOnly,
			SpeclibToParquetConverter converter) throws Exception {
		ProgressCallback cb = callback == null ? new ProgressCallback.NoOp() : callback;
		if (magDb == null) {
			if (logger != null) logger.error("Marker4 metadata: magDb is null");
			return false;
		}
		if (diaPar == null) {
			if (logger != null) logger.error("Marker4 metadata: DIA parameters are null");
			return false;
		}
		if (converter == null) {
			if (logger != null) logger.error("Marker4 metadata: speclib converter is null");
			return false;
		}

		File speclib = magDb.getMarker4LibFile();
		if (speclib == null || !speclib.exists()) {
			if (logger != null) logger.error("Marker4 metadata: canonical speclib missing");
			return false;
		}
		File libFolder = magDb.getMarker4LibraryFolder();
		if (libFolder == null) {
			if (logger != null) logger.error("Marker4 metadata: library folder null");
			return false;
		}
		if (!magDb.ensureMarker4LibraryFolderExists()) {
			if (logger != null) logger.error("Marker4 metadata: could not create library folder {}", libFolder.getAbsolutePath());
			return false;
		}

		String name = MagDbItem.MARKER4;
		File parquet = new File(libFolder, name + ".parquet");
		File infoFile = new File(libFolder, name + "_info.txt");

		DiannParameter diannPar = buildDiannParForLibMetadata(diaPar, speclibOnly);
		if (!parquet.exists()) {
			cb.onMessage("Marker4: converting speclib to parquet (" + parquet.getName() + ")");
			converter.convert(diannPar, speclib, parquet, libFolder, cb, logger);
		}
		if (!parquet.exists()) {
			if (logger != null) logger.error("Marker4 metadata: parquet conversion failed, expected {}", parquet.getAbsolutePath());
			return false;
		}

		HashSet<String> pepSet = DiaNNParquetReader.getPeptideSet(parquet.getAbsolutePath());
		int n = pepSet != null ? pepSet.size() : 0;
		try (PrintWriter writer = new PrintWriter(infoFile)) {
			writer.println("Data set size:\t" + n);
			writer.println("Enzyme:\t" + ENZYME);
			writer.println("Miss cleavages:\t" + MISS);
			writer.println("MatchToMag:\t" + MATCH_TO_MAG);
		}
		magDb.addPepLib(name, n, ENZYME, MISS, MATCH_TO_MAG);
		return true;
	}

	public static boolean ensureMarker4LibraryMetadata(File diannExe, File condaPrefix, MagDbItem magDb,
			MetaParameterDia diaPar, ProgressCallback callback, Logger logger) throws Exception {
		return ensureMarker4LibraryMetadata(magDb, diaPar, callback, logger, true,
				(diannPar, speclib, parquet, libFolder, cb, log) -> {
					String cmd = diannPar.generateSpeclib2Parquet(speclib.getAbsolutePath(),
							parquet.getAbsolutePath());
					int exit = DiaNNHeadlessRunner.run(diannExe, cmd, libFolder, cb, condaPrefix);
					if (exit != 0 && log != null) {
						log.error("Marker4 metadata: DIA-NN speclib->parquet exit code {}", exit);
					}
					return exit == 0;
				});
	}
}
