package uk.ac.quadram.metapilot.core.runner;

import java.io.File;
import java.io.IOException;
import java.io.InterruptedIOException;
import java.io.RandomAccessFile;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;

import uk.ac.quadram.metapilot.core.workflow.ProgressCallback;
import uk.ac.quadram.metapilot.spectra.io.MSConvertor;
import uk.ac.quadram.metapilot.task.dda.mag.par.DDAMagOpenSearchParameters;
import uk.ac.quadram.metapilot.task.par.MetaData;

/**
 * Converts Thermo {@code .raw} input to mzML before the open search starts;
 * mzML input is left exactly as loaded.
 * <p>
 * Crystal-C, stage-1 MSBooster, Sage and the variant MS2 ranking read only
 * mzML, so every step reads the converted file. The mzML go to the result
 * folder's {@code spectra/} ({@code MetaBaseParameter.getSpectraFile()}), where
 * the other workflows put converted spectra. The workflow reads the spectra
 * from one folder, so a mix of {@code .raw} and other input is refused.
 * <p>
 * msconvert's own output goes to the progress callback, which writes the run's
 * {@code logging.log}. A failed file is reported there with msconvert's
 * message and the ProteoWizard release, and the run stops once every other
 * file has been converted, so a rerun converts only what failed.
 */
final class RawSpectraConversion {

	/** msconvert's wording for an instrument its release does not know. */
	private static final String[] VERSION_PROBLEM_MARKERS = { "unable to parse instrument model",
			"latest version of ProteoWizard" };

	/** Last element msconvert writes; without it the file was cut short. */
	private static final String MZML_END_TAG = "</indexedmzML>";

	private RawSpectraConversion() {
	}

	/**
	 * Returns the experiment rows the workflow should read. With no
	 * {@code .raw} input this is {@code params.getExperimentData()}, untouched.
	 * Otherwise each raw file is converted into {@code params.getSpectraFile()}
	 * unless a complete conversion no older than it is already there, and both
	 * the experiment rows and the metadata FlashLFQ reads are pointed at the
	 * mzML.
	 */
	static List<String[]> prepare(DDAMagOpenSearchParameters params, ProgressCallback cb) throws IOException {
		List<String[]> experiments = params.getExperimentData();
		List<String> others = new ArrayList<>();
		boolean anyRaw = false;
		for (String[] experiment : experiments) {
			if (isRaw(experiment[0])) {
				anyRaw = true;
			} else {
				others.add(experiment[0]);
			}
		}
		if (!anyRaw) {
			return experiments;
		}
		if (!others.isEmpty()) {
			throw new IOException("The input mixes Thermo .raw files with other spectra (" + others.get(0)
					+ (others.size() > 1 ? " and " + (others.size() - 1) + " more" : "")
					+ "). The workflow reads all spectra from one folder, so load either only .raw files, "
					+ "which are converted to mzML, or only mzML files.");
		}

		File outDir = params.getSpectraFile();
		Map<String, String> samples = new HashMap<>();
		Map<String, String> converted = new HashMap<>();
		List<File> pending = new ArrayList<>();
		for (String[] experiment : experiments) {
			File raw = new File(experiment[0]);
			String sample = baseName(raw);
			String clash = samples.put(sample.toLowerCase(Locale.ROOT), experiment[0]);
			if (clash != null) {
				throw new IOException("Two raw files share the sample name " + sample + " (" + clash + ", "
						+ experiment[0] + "); their mzML would overwrite each other.");
			}
			File mzml = mzmlFor(raw, outDir);
			converted.put(experiment[0], mzml.getAbsolutePath());
			if (isCurrent(mzml, raw)) {
				cb.onMessage("Spectra: " + mzml.getName() + " already converted from " + raw.getName() + ", reusing");
			} else if (!raw.isFile()) {
				throw new IOException("Raw file not found: " + raw.getAbsolutePath());
			} else {
				pending.add(raw);
			}
		}
		if (!pending.isEmpty()) {
			convert(params, pending, outDir, cb);
		}

		List<String[]> rows = new ArrayList<>(experiments.size());
		for (String[] experiment : experiments) {
			String[] row = experiment.clone();
			row[0] = converted.get(experiment[0]);
			rows.add(row);
		}
		params.setExperimentData(rows);
		MetaData metadata = params.getMetadata();
		if (metadata != null && metadata.getRawFiles() != null) {
			String[] files = metadata.getRawFiles().clone();
			for (int i = 0; i < files.length; i++) {
				files[i] = converted.getOrDefault(files[i], files[i]);
			}
			metadata.setRawFiles(files);
		}
		cb.onMessage("Spectra: " + rows.size() + " Thermo raw file(s) are read as mzML from "
				+ outDir.getAbsolutePath());
		return rows;
	}

	private static void convert(DDAMagOpenSearchParameters params, List<File> pending, File outDir,
			ProgressCallback cb) throws IOException {
		String exe = params.getMsconvertExecutable();
		if (exe == null || exe.isBlank() || !new File(exe).isFile()) {
			throw new IOException(pending.size() + " Thermo .raw file(s) must be converted to mzML first, but msconvert "
					+ (exe == null || exe.isBlank() ? "is not configured" : "was not found at " + exe)
					+ ". Set Resources > msconvert to ProteoWizard's msconvert.exe, or load mzML files.");
		}
		String version = MSConvertor.version(exe);
		String release = version == null ? "ProteoWizard release unknown" : version;
		cb.onMessage("Spectra: converting " + pending.size() + " Thermo raw file(s) to mzML with " + exe + " ("
				+ release + "): " + String.join(" ", MSConvertor.raw2mzMLCommand(exe, pending.get(0), outDir)));

		int threads = Math.min(pending.size(), Math.max(1, params.getThreadCount()));
		ExecutorService pool = Executors.newFixedThreadPool(threads);
		List<String> failures = new ArrayList<>();
		try {
			List<Future<String>> results = new ArrayList<>(pending.size());
			for (File raw : pending) {
				results.add(pool.submit(() -> convertOne(exe, release, raw, outDir, cb)));
			}
			for (int i = 0; i < pending.size(); i++) {
				String problem;
				try {
					problem = results.get(i).get();
				} catch (ExecutionException e) {
					Throwable cause = e.getCause() == null ? e : e.getCause();
					problem = cause.getClass().getSimpleName() + ": " + cause.getMessage();
				}
				if (problem != null) {
					failures.add(pending.get(i).getName() + " - " + problem);
				}
			}
		} catch (InterruptedException e) {
			Thread.currentThread().interrupt();
			throw new InterruptedIOException("Interrupted while converting Thermo raw files to mzML");
		} finally {
			pool.shutdownNow();
		}

		if (!failures.isEmpty()) {
			for (String failure : failures) {
				cb.onMessage("ERROR: msconvert did not convert " + failure);
			}
			throw new IOException("msconvert did not convert " + failures.size() + " of " + pending.size()
					+ " Thermo raw file(s) to mzML (" + release + "); msconvert's messages are in logging.log. "
					+ "First: " + failures.get(0));
		}
	}

	/**
	 * Converts {@code raw} straight into {@code outDir}. A file left by a
	 * failed or interrupted conversion is deleted, so it is never mistaken for
	 * a finished one.
	 *
	 * @return null on success, otherwise why {@code raw} was not converted
	 */
	private static String convertOne(String exe, String release, File raw, File outDir, ProgressCallback cb)
			throws IOException, InterruptedException {
		String sample = baseName(raw);
		File target = mzmlFor(raw, outDir);
		Files.deleteIfExists(target.toPath());
		List<String> output = new ArrayList<>();
		long start = System.nanoTime();
		boolean done = false;
		try {
			int exit = MSConvertor.raw2mzML(exe, raw, outDir, line -> {
				output.add(line);
				cb.onMessage("msconvert [" + sample + "]: " + line);
			});
			if (exit != 0 || !isComplete(target)) {
				String problem = "exit code " + exit + (isComplete(target) ? "" : ", no complete mzML written");
				if (isVersionProblem(output)) {
					cb.onMessage("WARNING: " + raw.getName() + " was acquired on an instrument this msconvert ("
							+ release + ") does not recognise. This is a ProteoWizard version problem: install a "
							+ "newer ProteoWizard and point Resources > msconvert at its msconvert.exe.");
					problem += "; ProteoWizard version problem (" + release + " does not recognise the instrument)";
				}
				return problem;
			}
			done = true;
			cb.onMessage("Spectra: converted " + raw.getName() + " to " + target.getName() + " in "
					+ Math.round((System.nanoTime() - start) / 1e9) + " s");
			return null;
		} finally {
			if (!done) {
				Files.deleteIfExists(target.toPath());
			}
		}
	}

	/** A complete conversion that is no older than its raw file. */
	static boolean isCurrent(File mzml, File raw) {
		return isComplete(mzml) && (!raw.exists() || mzml.lastModified() >= raw.lastModified());
	}

	/** True when {@code mzml} ends with msconvert's closing tag. */
	static boolean isComplete(File mzml) {
		if (!mzml.isFile() || mzml.length() == 0) {
			return false;
		}
		try (RandomAccessFile in = new RandomAccessFile(mzml, "r")) {
			int n = (int) Math.min(256, in.length());
			byte[] tail = new byte[n];
			in.seek(in.length() - n);
			in.readFully(tail);
			return new String(tail, StandardCharsets.US_ASCII).contains(MZML_END_TAG);
		} catch (IOException e) {
			return false;
		}
	}

	static boolean isVersionProblem(List<String> msconvertOutput) {
		for (String line : msconvertOutput) {
			for (String marker : VERSION_PROBLEM_MARKERS) {
				if (line.contains(marker)) {
					return true;
				}
			}
		}
		return false;
	}

	static boolean isRaw(String path) {
		return path != null && path.toLowerCase(Locale.ROOT).endsWith(".raw");
	}

	static File mzmlFor(File raw, File outDir) {
		return new File(outDir, baseName(raw) + ".mzML");
	}

	private static String baseName(File file) {
		String name = file.getName();
		int dot = name.lastIndexOf('.');
		return dot > 0 ? name.substring(0, dot) : name;
	}
}
