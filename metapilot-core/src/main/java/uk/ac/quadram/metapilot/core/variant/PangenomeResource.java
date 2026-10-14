package uk.ac.quadram.metapilot.core.variant;

import java.io.File;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.util.function.Consumer;

import org.json.JSONObject;

/**
 * Where the prepared UHGP-100 files are, and whether they are complete.
 * <p>
 * Resolution, first match wins:
 * <ol>
 * <li>{@code uhgp100_path} set to a directory: the prepared files in it, or a
 * {@code uhgp-100.tar.gz} in it to prepare from.</li>
 * <li>{@code uhgp100_path} set to {@code uhgp-100.faa.gz}: that file, with the
 * species table and manifest beside it.</li>
 * <li>{@code uhgp100_path} set to {@code uhgp-100.tar.gz}: prepared once into
 * its own directory.</li>
 * <li>Unset: {@code protein_catalogue/} under the selected MAG database
 * version, which is where MGnify publishes the archive relative to
 * {@code genomes-all_metadata.tsv}.</li>
 * </ol>
 * The manifest written by the prepare step records the output sizes. A
 * {@code uhgp-100.faa.gz} whose size disagrees with it is an interrupted or
 * replaced preparation and is not used.
 */
public final class PangenomeResource {

	public static final String FAA = "uhgp-100.faa.gz";
	public static final String REP_SPECIES = "uhgp-100.rep_species.tsv.gz";
	public static final String MANIFEST = "uhgp-100.prepare_manifest.json";
	public static final String ARCHIVE = "uhgp-100.tar.gz";
	public static final String METADATA = "genomes-all_metadata.tsv";
	public static final String DEFAULT_SUBDIR = "protein_catalogue";

	public final File directory;
	public final File faaGz;
	/** Null when absent: hits then count only the representative's own species. */
	public final File repSpeciesGz;
	/** Null when absent: species scope is unavailable and only ALL is reported. */
	public final File metadataTsv;
	/** Non-null when the FASTA still has to be prepared from this archive. */
	public final File archiveToPrepare;

	private PangenomeResource(File directory, File faaGz, File repSpeciesGz, File metadataTsv,
			File archiveToPrepare) {
		this.directory = directory;
		this.faaGz = faaGz;
		this.repSpeciesGz = repSpeciesGz;
		this.metadataTsv = metadataTsv;
		this.archiveToPrepare = archiveToPrepare;
	}

	/**
	 * Locates the resource without touching it; null when there is none or it is
	 * incomplete with nothing to prepare it from.
	 *
	 * @param magVersionDir the MAG database version directory
	 *                      ({@code MagDbItem.getCurrentFile()}), or null
	 */
	public static PangenomeResource locate(String configuredPath, File magVersionDir) {
		File dir = null;
		File faa = null;
		File archive = null;
		if (configuredPath != null && !configuredPath.isBlank()) {
			File f = new File(configuredPath.trim());
			if (f.isDirectory()) {
				dir = f;
			} else if (f.isFile()) {
				String name = f.getName().toLowerCase(java.util.Locale.ROOT);
				dir = f.getAbsoluteFile().getParentFile();
				if (name.endsWith(".tar.gz") || name.endsWith(".tgz")) {
					archive = f;
				} else {
					faa = f;
				}
			} else {
				return null;
			}
		} else if (magVersionDir != null) {
			File d = new File(magVersionDir, DEFAULT_SUBDIR);
			if (!d.isDirectory()) {
				return null;
			}
			dir = d;
		} else {
			return null;
		}
		if (faa == null) {
			faa = new File(dir, FAA);
		}
		if (archive == null) {
			File a = new File(dir, ARCHIVE);
			archive = a.isFile() ? a : null;
		}
		File metadata = firstFile(new File(dir, METADATA), dir.getParentFile() == null ? null
				: new File(dir.getParentFile(), METADATA),
				magVersionDir == null ? null : new File(magVersionDir, METADATA));
		boolean complete = faa.isFile() && matchesManifest(faa, new File(dir, MANIFEST));
		if (!complete && archive == null) {
			return null;
		}
		File rep = new File(dir, REP_SPECIES);
		return new PangenomeResource(dir, faa, complete && rep.isFile() ? rep : null, metadata,
				complete ? null : archive);
	}

	public boolean needsPrepare() {
		return archiveToPrepare != null;
	}

	/** Runs the one-time prepare step and returns the resulting resource. */
	public PangenomeResource prepare(Consumer<String> progress) throws IOException {
		if (archiveToPrepare == null) {
			return this;
		}
		PangenomeResourcePreparer.prepare(archiveToPrepare, metadataTsv, directory, progress);
		File rep = new File(directory, REP_SPECIES);
		return new PangenomeResource(directory, new File(directory, FAA), rep.isFile() ? rep : null, metadataTsv,
				null);
	}

	public String describe() {
		return faaGz.getAbsolutePath() + (repSpeciesGz == null ? " (no species table)" : "")
				+ (metadataTsv == null ? " (no genome metadata)" : "");
	}

	/** True when there is no manifest (hand-placed file) or its recorded size matches. */
	private static boolean matchesManifest(File faa, File manifest) {
		if (!manifest.isFile()) {
			return true;
		}
		try {
			JSONObject m = new JSONObject(new String(Files.readAllBytes(manifest.toPath()), StandardCharsets.UTF_8));
			long expected = m.optLong("out_faa_gz_bytes", -1L);
			return expected < 0 || expected == faa.length();
		} catch (IOException | RuntimeException e) {
			return false;
		}
	}

	private static File firstFile(File... candidates) {
		for (File f : candidates) {
			if (f != null && f.isFile()) {
				return f;
			}
		}
		return null;
	}
}
