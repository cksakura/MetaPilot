package uk.ac.quadram.metapilot.core.runner;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import uk.ac.quadram.metapilot.core.dia.PeptideListLibrarySupport;
import uk.ac.quadram.metapilot.core.diann.DiaNNHeadlessRunner;
import uk.ac.quadram.metapilot.core.workflow.ProgressCallback;
import uk.ac.quadram.metapilot.dbSearch.diann.DiaNNParquetReader;
import uk.ac.quadram.metapilot.dbSearch.diann.DiannParameter;
import uk.ac.quadram.metapilot.task.dia.par.MetaParameterDia;
import uk.ac.quadram.metapilot.task.dia.par.MetaParameterDiaMarker;
import uk.ac.quadram.metapilot.task.mag.MagDbItem;

/**
 * Headless DIA library build into the MAG database ({@code libraries/<name>/}), the steps of the GUI library task
 * ({@code DiaLibCreateTask}) for a peptide list or a protein FASTA:
 * <ol>
 *   <li>peptide list: the unique unmodified peptides are written to {@code <name>.fasta}, one sorted entry each, and
 *   the peptide-entry rules of {@link PeptideListLibrarySupport} apply (no protein N-terminal acetylation or Met
 *   excision, missed cleavages that keep the peptides intact); protein FASTA: used as it is;</li>
 *   <li>DIA-NN predicts {@code <name>.speclib} (enzyme Trypsin/P, the modification settings of the parameter file);</li>
 *   <li>the library is converted to {@code <name>.parquet};</li>
 *   <li>{@code <name>_info.txt} records size, enzyme, missed cleavages, Microbiome/Standalone type and, for a peptide
 *   list, the peptide-entry details; {@code dia-seed} finds the library through this file.</li>
 * </ol>
 */
public final class DiaLibraryHeadlessRunner {

	private static final Logger LOGGER = LogManager.getLogger(DiaLibraryHeadlessRunner.class);
	private static final String TASK_NAME = "DIA library";
	private static final String ENZYME = "Trypsin/P";

	public enum Source { PEPTIDES, PROTEINS }

	private DiaLibraryHeadlessRunner() {
	}

	public static final class Result {
		public final File folder;
		public final File speclib;
		public final File parquet;
		public final File infoFile;
		public final int dataSetSize;
		public final int missedCleavages;

		Result(File folder, File speclib, File parquet, File infoFile, int dataSetSize, int missedCleavages) {
			this.folder = folder;
			this.speclib = speclib;
			this.parquet = parquet;
			this.infoFile = infoFile;
			this.dataSetSize = dataSetSize;
			this.missedCleavages = missedCleavages;
		}
	}

	/**
	 * @param params                   DIA parameter file (Seed or Marker): MAG database, threads, modifications
	 * @param input                    peptide table (Base sequence / Sequence column) or peptide FASTA for
	 *                                 {@link Source#PEPTIDES}; protein FASTA for {@link Source#PROTEINS}
	 * @param requestedMissedCleavages missed cleavages for proteins; for peptides the minimum (raised to keep the
	 *                                 peptides intact)
	 * @param matchToMag               true = Microbiome library (peptides mapped to MAG proteins), false = Standalone
	 * @param overwrite                rebuild when the library folder already holds a library
	 */
	public static Result build(MetaParameterDia params, File magDbConfigJson, File diannExe, String libName, File input,
			Source source, int requestedMissedCleavages, boolean matchToMag, boolean overwrite, ProgressCallback cb,
			File condaPrefix) throws Exception {
		if (cb == null) cb = new ProgressCallback.NoOp();
		if (params == null) throw new IllegalArgumentException("params is null");
		if (diannExe == null || !diannExe.exists()) {
			throw new IOException("DIA-NN not found: " + (diannExe == null ? "null" : diannExe.getAbsolutePath()));
		}
		if (libName == null || !libName.matches("[A-Za-z0-9._-]+")) {
			throw new IOException("Library name must contain only letters, digits, '.', '_' or '-': " + libName);
		}
		if (input == null || !input.isFile()) {
			throw new IOException("Input file not found: " + (input == null ? "null" : input.getAbsolutePath()));
		}
		MagDbItem magDb = params instanceof MetaParameterDiaMarker
				? DiaMarkerHeadlessRunner.resolveMagDb((MetaParameterDiaMarker) params, magDbConfigJson)
				: params.getMagDbItem();
		if (magDb == null || magDb.getLibrariesRootFolder() == null) {
			throw new IOException("MAG database in the params file (mag_database) could not be resolved");
		}
		File folder = new File(magDb.getLibrariesRootFolder(), libName);
		File speclib = new File(folder, libName + ".speclib");
		File parquet = new File(folder, libName + ".parquet");
		File infoFile = new File(folder, libName + "_info.txt");
		File peptideFasta = new File(folder, libName + ".fasta");
		if (speclib.exists() || parquet.exists()) {
			if (!overwrite) {
				throw new IOException("Library " + libName + " already exists in " + folder.getAbsolutePath()
						+ " (use --overwrite to rebuild it)");
			}
			for (File f : new File[] { speclib, parquet, infoFile, peptideFasta,
					new File(folder, libName + ".predicted.speclib") }) {
				Files.deleteIfExists(f.toPath());
			}
		}
		if (!folder.isDirectory() && !folder.mkdirs()) {
			throw new IOException("Could not create library folder: " + folder.getAbsolutePath());
		}
		cb.onMessage(TASK_NAME + ": building " + libName + " in " + folder.getAbsolutePath());

		DiannParameter diannPar = new DiannParameter();
		diannPar.setThreads(params.getThreadCount());
		diannPar.setEnzyme(ENZYME);
		diannPar.setMissed_cleavages(requestedMissedCleavages);
		diannPar.setM_excision(params.isMExcision());
		diannPar.setCarbami_C(params.isCarbamiC());
		diannPar.setOxidation_M(params.isOxidationM());
		diannPar.setAcety_proN(params.isAcetyProN());
		diannPar.setPhospho_STY(params.isPhosphoSTY());
		diannPar.setK_GG(params.isKGG());

		File libFasta;
		int missed = requestedMissedCleavages;
		List<String> extraInfo = new ArrayList<>();
		if (source == Source.PEPTIDES) {
			String name = input.getName().toLowerCase();
			PeptideListLibrarySupport.PeptideList list = name.endsWith(".fasta") || name.endsWith(".fa")
					? PeptideListLibrarySupport.readPeptideFasta(input)
					: PeptideListLibrarySupport.readPeptideTable(input);
			if (list.peptides.isEmpty()) throw new IOException("No peptides read from " + input.getAbsolutePath());
			PeptideListLibrarySupport.writePeptideFasta(list, peptideFasta);
			missed = PeptideListLibrarySupport.configureForPeptideLibrary(diannPar, requestedMissedCleavages, list);
			extraInfo.addAll(PeptideListLibrarySupport.describe(list, diannPar));
			libFasta = peptideFasta;
		} else {
			libFasta = input;
			extraInfo.add("Library source:\tprotein FASTA");
		}
		extraInfo.add("Input file:\t" + input.getAbsolutePath());
		for (String line : extraInfo) cb.onMessage(TASK_NAME + ": " + line.replace('\t', ' '));

		cb.onMessage(TASK_NAME + ": predicting " + speclib.getName() + " from " + libFasta.getName());
		int exit = DiaNNHeadlessRunner.run(diannExe,
				diannPar.generateFasta(libFasta.getAbsolutePath(), speclib.getAbsolutePath(), false), folder, cb,
				condaPrefix);
		if (exit != 0) LOGGER.warn("{}: DIA-NN library prediction exit code {}", TASK_NAME, exit);
		File predicted = new File(folder, libName + ".predicted.speclib");
		if (!speclib.exists() && predicted.exists()) Files.move(predicted.toPath(), speclib.toPath());
		if (!speclib.exists()) throw new IOException("DIA-NN did not produce " + speclib.getAbsolutePath());

		cb.onMessage(TASK_NAME + ": converting to " + parquet.getName());
		exit = DiaNNHeadlessRunner.run(diannExe,
				diannPar.generateSpeclib2Parquet(speclib.getAbsolutePath(), parquet.getAbsolutePath()), folder, cb,
				condaPrefix);
		if (exit != 0) LOGGER.warn("{}: DIA-NN parquet conversion exit code {}", TASK_NAME, exit);
		if (!parquet.exists()) throw new IOException("DIA-NN did not produce " + parquet.getAbsolutePath());

		HashSet<String> pepSet = DiaNNParquetReader.getPeptideSet(parquet.getAbsolutePath());
		try (PrintWriter w = new PrintWriter(infoFile)) {
			w.println("Data set size:\t" + pepSet.size());
			w.println("Enzyme:\t" + ENZYME);
			w.println("Miss cleavages:\t" + missed);
			w.println("MatchToMag:\t" + matchToMag);
			for (String line : extraInfo) w.println(line);
			w.println("Built by:\tmetapilot-cli dia-library");
		}
		magDb.addPepLib(libName, pepSet.size(), ENZYME, missed, matchToMag);
		cb.onComplete("DIA library complete: " + folder.getAbsolutePath() + " (" + pepSet.size() + " peptides, "
				+ (matchToMag ? "Microbiome" : "Standalone") + ")");
		return new Result(folder, speclib, parquet, infoFile, pepSet.size(), missed);
	}
}
