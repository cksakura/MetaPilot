package uk.ac.quadram.metapilot.task.dia;

import java.io.File;
import java.io.IOException;
import java.util.Date;

import uk.ac.quadram.metapilot.core.runner.DiaMagPostAnalysisTsvExporter;
import uk.ac.quadram.metapilot.core.workflow.ProgressCallback;
import uk.ac.quadram.metapilot.task.mag.MagDbItem;

/**
 * Runnable Eclipse test/helper for DIA Seed post-analysis when the combined
 * DIA-NN result has not yet been mapped to MAG proteins.
 *
 * Input:
 * - combined result folder containing combined.parquet and combined.pr_matrix.tsv
 *   (or a parent folder containing combined/combined.parquet and
 *   combined/combined.pr_matrix.tsv)
 * - MAG DB current/version folder containing predicted/*.predicted.tsv and the
 *   taxonomy file
 *
 * Output:
 * - final_peptides.tsv
 * - razor_proteins.tsv
 * - all_proteins.tsv
 * - taxonomy_analysis/Genome.tsv
 * - taxonomy_analysis/Taxa_report.tsv
 */
public class DiaSeedMappedCombinedResultTest {

	private static final String COMBINED_SEARCH_RESULT_FOLDER = "Z:\\Test\\MetaPilot\\target_library\\combined";
	private static final String MAG_DB_CURRENT_FOLDER = "C:\\Users\\daz25hak\\MetaProAI\\database\\MAG\\human-gut\\v2.0.2";
	private static final String OUTPUT_RESULT_FOLDER = "Z:\\Test\\MetaPilot\\target_library\\result";

	private static final String MAG_CATALOGUE_ID = "human-gut";
	private static final String MAG_IDENTIFIER = "MGYG";
	private static final String ORIGINAL_DB_FOLDER_NAME = "original_db";
	private static final String FUNCTION_FOLDER_NAME = "eggNOG";
	private static final String TAXA_FILE_NAME = "genomes-all_metadata.tsv";
	private static final String GTDB_NCBI_MAPPING_TSV = "";
	private static final String FUNC_DEF_DB = "C:\\Users\\daz25hak\\MetaPilot\\MetaPilot\\Resources\\function\\func_def.db";

	public static void main(String[] args) throws Exception {
		File input = new File(value(args, 0, COMBINED_SEARCH_RESULT_FOLDER));
		File magDbCurrent = new File(value(args, 1, MAG_DB_CURRENT_FOLDER));
		File resultDir = resolveResultDir(input, value(args, 2, OUTPUT_RESULT_FOLDER));
		File combinedDir = resolveCombinedDir(input);

		File combinedParquet = new File(combinedDir, "combined.parquet");
		File combinedPrMatrix = new File(combinedDir, "combined.pr_matrix.tsv");
		requireFile(combinedParquet, "combined parquet");
		requireFile(combinedPrMatrix, "combined peptide matrix");
		requireDirectory(new File(magDbCurrent, "predicted"), "MAG predicted peptide folder");

		MagDbItem magDb = buildMagDbItem(magDbCurrent);
		File gtdbNcbiMappingTsv = optionalFile(value(args, 3, GTDB_NCBI_MAPPING_TSV));
		File funcDefDb = optionalFile(value(args, 4, FUNC_DEF_DB));

		System.out.println("DIA Seed mapped post-analysis");
		System.out.println("  combined folder: " + combinedDir.getAbsolutePath());
		System.out.println("  MAG DB folder:   " + magDbCurrent.getAbsolutePath());
		System.out.println("  output folder:   " + resultDir.getAbsolutePath());

		DiaMagPostAnalysisTsvExporter.Result result = DiaMagPostAnalysisTsvExporter.export(
				resultDir,
				magDb,
				combinedParquet,
				combinedPrMatrix,
				gtdbNcbiMappingTsv,
				"SEED",
				funcDefDb,
				MAG_CATALOGUE_ID,
				"Intensity",
				new ConsoleProgress());

		System.out.println("Finished.");
		System.out.println("  peptides: " + result.finalPeptidesTsv.getAbsolutePath());
		System.out.println("  proteins: " + result.razorProteinsTsv.getAbsolutePath());
		System.out.println("  genomes:  " + result.genomeTsv.getAbsolutePath());
	}

	private static String value(String[] args, int index, String fallback) {
		return args != null && args.length > index && args[index] != null && !args[index].isBlank()
				? args[index]
				: fallback;
	}

	private static File resolveCombinedDir(File input) throws IOException {
		if (new File(input, "combined.parquet").isFile() && new File(input, "combined.pr_matrix.tsv").isFile()) {
			return input.getCanonicalFile();
		}
		File combined = new File(input, "combined");
		if (new File(combined, "combined.parquet").isFile() && new File(combined, "combined.pr_matrix.tsv").isFile()) {
			return combined.getCanonicalFile();
		}
		throw new IOException("Cannot find combined.parquet and combined.pr_matrix.tsv under: "
				+ input.getAbsolutePath());
	}

	private static File resolveResultDir(File input, String configuredOutput) throws IOException {
		if (configuredOutput != null && !configuredOutput.isBlank()) {
			File out = new File(configuredOutput);
			if (!out.exists() && !out.mkdirs()) {
				throw new IOException("Cannot create output folder: " + out.getAbsolutePath());
			}
			return out.getCanonicalFile();
		}
		File combinedDir = resolveCombinedDir(input);
		File parent = combinedDir.getParentFile();
		return (parent == null ? combinedDir : parent).getCanonicalFile();
	}

	private static MagDbItem buildMagDbItem(File currentFolder) throws IOException {
		requireDirectory(currentFolder, "MAG DB current folder");
		File repository = currentFolder.getParentFile();
		if (repository == null) {
			repository = currentFolder;
		}
		MagDbItem magDb = new MagDbItem(MAG_CATALOGUE_ID, MAG_IDENTIFIER, repository);
		magDb.addVersion(currentFolder.getName(), 0, new Date(), true, ORIGINAL_DB_FOLDER_NAME, FUNCTION_FOLDER_NAME,
				TAXA_FILE_NAME);
		magDb.setUsedVersion(magDb.getVersionList().get(0));
		magDb.setCurrentFile(currentFolder);
		return magDb;
	}

	private static File optionalFile(String path) {
		if (path == null || path.isBlank()) {
			return null;
		}
		File file = new File(path);
		return file.isFile() ? file : null;
	}

	private static void requireFile(File file, String label) throws IOException {
		if (file == null || !file.isFile()) {
			throw new IOException(label + " not found: " + (file == null ? "null" : file.getAbsolutePath()));
		}
	}

	private static void requireDirectory(File file, String label) throws IOException {
		if (file == null || !file.isDirectory()) {
			throw new IOException(label + " not found: " + (file == null ? "null" : file.getAbsolutePath()));
		}
	}

	private static final class ConsoleProgress implements ProgressCallback {
		@Override
		public void onProgress(double workDone, double max) {
			if (max > 0) {
				System.out.println(String.format("Progress: %.1f%%", workDone * 100.0 / max));
			}
		}

		@Override
		public void onMessage(String message) {
			System.out.println(message);
		}

		@Override
		public void onComplete(Object result) {
			System.out.println("Complete: " + result);
		}

		@Override
		public void onError(Exception error) {
			error.printStackTrace();
		}
	}
}
