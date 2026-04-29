package uk.ac.quadram.metapilot.task.dia;

import java.io.File;
import java.io.PrintWriter;
import java.util.HashSet;

import org.apache.logging.log4j.Logger;

import uk.ac.quadram.metapilot.dbSearch.diann.DiaNNParquetReader;
import uk.ac.quadram.metapilot.dbSearch.diann.DiaNNTask;
import uk.ac.quadram.metapilot.dbSearch.diann.DiannParameter;
import uk.ac.quadram.metapilot.resources.MetaResourceParameters;
import uk.ac.quadram.metapilot.resources.MetaResourceParameters.ToolType;
import uk.ac.quadram.metapilot.task.dia.par.MetaParameterDia;
import uk.ac.quadram.metapilot.task.mag.MagDbItem;

/**
 * Marker4 spectral library layout under {@code libraries/Marker4/}: matches {@link DiaLibCreateTask} (parquet +
 * {@code Marker4_info.txt} + {@link MagDbItem#addPepLib}).
 */
public final class Marker4LibrarySupport {

	private static final String ENZYME = "Trypsin/P";
	private static final int MISS = 1;
	private static final boolean MATCH_TO_MAG = true;

	private Marker4LibrarySupport() {
	}

	private static DiannParameter buildDiannParForLibMetadata(MetaParameterDia diaPar) {
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
		return diannPar;
	}

	/**
	 * Ensures {@code Marker4.parquet}, {@code Marker4_info.txt}, and in-memory registration exist for the canonical
	 * speclib (same steps as the end of {@link DiaLibCreateTask}).
	 */
	public static boolean ensureMarker4LibraryMetadata(MagDbItem magDb, MetaParameterDia diaPar,
			MetaResourceParameters mrp, Logger logger) {
		File speclib = magDb.getMarker4LibFile();
		if (speclib == null || !speclib.exists()) {
			logger.error("Marker4 metadata: canonical speclib missing");
			return false;
		}
		File libFolder = magDb.getMarker4LibraryFolder();
		if (libFolder == null) {
			logger.error("Marker4 metadata: library folder null");
			return false;
		}
		if (!magDb.ensureMarker4LibraryFolderExists()) {
			logger.error("Marker4 metadata: could not create library folder {}", libFolder.getAbsolutePath());
			return false;
		}
		String name = MagDbItem.MARKER4;
		File parquet = new File(libFolder, name + ".parquet");
		File infoFile = new File(libFolder, name + "_info.txt");

		DiannParameter diannPar = buildDiannParForLibMetadata(diaPar);
		DiaNNTask diaNNTask = new DiaNNTask(mrp.getToolPath(ToolType.DIANN));

		if (!parquet.exists()) {
			diaNNTask.addTask(diannPar, speclib.getAbsolutePath(), parquet.getAbsolutePath());
			diaNNTask.run(1);
		}
		if (!parquet.exists()) {
			logger.error("Marker4 metadata: parquet conversion failed, expected {}", parquet.getAbsolutePath());
			return false;
		}

		try {
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
		} catch (Exception e) {
			logger.error("Marker4 metadata: failed writing info or registering library", e);
			return false;
		}
	}
}
