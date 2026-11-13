package uk.ac.quadram.metapilot.core.opensearch;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.StandardCopyOption;
import java.util.ArrayList;
import java.util.List;
import java.util.Locale;

/**
 * Which file the open-search workflow reads spectra from, for each input run.
 * <p>
 * An mzML input is read directly. A vendor input (Thermo {@code .raw}, Bruker {@code .d}) is searched
 * by MSFragger as it is, and read through the {@code <run>_uncalibrated.mzML} MSFragger writes during
 * the marker search ({@code write_uncalibrated_mzml}). MSFragger's scan numbers refer to that file,
 * and for timsTOF ddaPASEF it holds one MS2 spectrum per precursor, summed over its PASEF frames,
 * which a generic converter does not. The file is moved to {@code <resultRoot>/spectra_mzml/<run>.mzML}
 * so every reader - MSBooster, Crystal-C, the MS2 features, shift localisation, FlashLFQ and the
 * spectra summary - finds it under the run's own name, the name the PSMs carry.
 * <p>
 * A run's inputs are either all mzML or all vendor files: the downstream tools take one spectra
 * folder.
 */
public final class OpenSearchSpectra {

	/** Under the engine result root: the mzML read for each vendor-format input. */
	public static final String CONVERTED_DIR = "spectra_mzml";

	/** MSFragger's name for the uncalibrated mzML it writes beside a vendor input. */
	static final String UNCALIBRATED_SUFFIX = "_uncalibrated.mzML";

	private OpenSearchSpectra() {
	}

	public static boolean isMzml(File input) {
		return input != null && input.getName().toLowerCase(Locale.ROOT).endsWith(".mzml");
	}

	/** True when the inputs are vendor files, which MSFragger must convert; false for mzML inputs. */
	public static boolean needsConversion(List<File> inputs) throws IOException {
		int mzml = 0;
		for (File input : inputs) {
			if (isMzml(input)) {
				mzml++;
			}
		}
		if (mzml > 0 && mzml < inputs.size()) {
			throw new IOException("The open search needs all inputs as mzML, or all as vendor files (.raw/.d); "
					+ "this run mixes " + mzml + " mzML with " + (inputs.size() - mzml) + " vendor files.");
		}
		return mzml == 0 && !inputs.isEmpty();
	}

	/** Run name: the file name without its extension ({@code X.mzML}, {@code X.raw}, {@code X.d} give X). */
	public static String runName(File input) {
		String name = input.getName();
		int dot = name.lastIndexOf('.');
		return dot > 0 ? name.substring(0, dot) : name;
	}

	public static File convertedDir(File resultRoot) {
		return new File(resultRoot, CONVERTED_DIR);
	}

	public static File convertedMzml(File resultRoot, File input) {
		return new File(convertedDir(resultRoot), runName(input) + ".mzML");
	}

	/** The mzML the workflow reads for one input: the input itself, or its converted copy. */
	public static File readable(File input, File resultRoot) {
		return isMzml(input) ? input : convertedMzml(resultRoot, input);
	}

	public static List<File> readable(List<File> inputs, File resultRoot) {
		List<File> out = new ArrayList<>(inputs.size());
		for (File input : inputs) {
			out.add(readable(input, resultRoot));
		}
		return out;
	}

	/** The folder the spectrum-reading tools are pointed at: the inputs' folder, or {@link #CONVERTED_DIR}. */
	public static File spectraDir(List<File> inputs, File resultRoot) throws IOException {
		if (inputs.isEmpty()) {
			return null;
		}
		return needsConversion(inputs) ? convertedDir(resultRoot) : inputs.get(0).getParentFile();
	}

	/** True when the input is mzML, or its converted mzML is already in place. */
	public static boolean isReadable(File input, File resultRoot) {
		File f = readable(input, resultRoot);
		return f.isFile() && f.length() > 0;
	}

	/** Where MSFragger writes the uncalibrated mzML for a vendor input: beside it. */
	static File uncalibratedOutput(File input) {
		return new File(input.getAbsoluteFile().getParentFile(), runName(input) + UNCALIBRATED_SUFFIX);
	}

	/**
	 * Moves the uncalibrated mzML MSFragger just wrote for a vendor input to its converted path,
	 * replacing an older copy. Throws when MSFragger wrote none.
	 */
	public static File collect(File input, File resultRoot) throws IOException {
		File produced = uncalibratedOutput(input);
		if (!produced.isFile() || produced.length() == 0) {
			throw new IOException("MSFragger wrote no uncalibrated mzML for " + input.getName() + " (expected "
					+ produced.getAbsolutePath() + "); the workflow reads the spectra of .raw and .d inputs from "
					+ "it");
		}
		File target = convertedMzml(resultRoot, input);
		target.getParentFile().mkdirs();
		// Written beside the final name and renamed last, so a crash never leaves a partial
		// mzML that a resumed run would take as complete.
		File partial = new File(target.getParentFile(), target.getName() + ".part");
		try {
			Files.move(produced.toPath(), partial.toPath(), StandardCopyOption.REPLACE_EXISTING);
		} catch (IOException crossDevice) {
			Files.copy(produced.toPath(), partial.toPath(), StandardCopyOption.REPLACE_EXISTING);
			Files.delete(produced.toPath());
		}
		Files.move(partial.toPath(), target.toPath(), StandardCopyOption.REPLACE_EXISTING);
		return target;
	}
}
