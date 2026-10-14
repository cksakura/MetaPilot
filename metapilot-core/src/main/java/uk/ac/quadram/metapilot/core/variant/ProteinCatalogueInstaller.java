package uk.ac.quadram.metapilot.core.variant;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.StandardCopyOption;
import java.util.Locale;

import uk.ac.quadram.metapilot.core.tools.ResumableDownload;
import uk.ac.quadram.metapilot.core.workflow.ProgressCallback;

/**
 * Downloads and prepares the UHGP-100 protein catalogue beside a MAG database,
 * as an optional part of the MAG database download.
 * <p>
 * Everything lands in {@code <catalogue>/<version>/protein_catalogue/}, which is
 * where {@link PangenomeResource#locate} looks when {@code uhgp100_path} is
 * unset, so a run finds it without configuration. The archive is downloaded
 * under a {@code .download} name and renamed only when complete: a partial
 * {@code uhgp-100.tar.gz} would otherwise look like something to prepare from.
 * <p>
 * Only the human-gut catalogue is supported; its archive name and layout were
 * verified for v2.0.2. Other MGnify catalogues publish protein catalogues under
 * other names.
 */
public final class ProteinCatalogueInstaller {

	public static final String SUPPORTED_CATALOGUE = "human-gut";

	/** {@code uhgp-100.tar.gz}, human-gut v2.0.2, Content-Length on 2026-10-07. */
	static final long ARCHIVE_BYTES = 28_342_719_422L;
	/** Prepared {@code uhgp-100.faa.gz} plus species table, measured 2026-10-07. */
	static final long PREPARED_BYTES = 22_872_894_230L;
	private static final String DOWNLOAD_SUFFIX = ".download";

	private ProteinCatalogueInstaller() {
	}

	public static boolean supports(String catalogueId) {
		return SUPPORTED_CATALOGUE.equalsIgnoreCase(catalogueId == null ? "" : catalogueId.trim());
	}

	public static File directory(File versionDir) {
		return new File(versionDir, PangenomeResource.DEFAULT_SUBDIR);
	}

	/** True when a complete prepared catalogue is already in place. */
	public static boolean isInstalled(File versionDir) {
		PangenomeResource r = PangenomeResource.locate(null, versionDir);
		return r != null && !r.needsPrepare() && r.faaGz.isFile();
	}

	/**
	 * Downloads (resuming a partial download) and prepares the catalogue. Does
	 * nothing but tidy up when it is already prepared.
	 *
	 * @param catalogueBaseUrl the catalogue version's URL, e.g.
	 *                         {@code http://ftp.ebi.ac.uk/pub/databases/metagenomics/mgnify_genomes/human-gut/v2.0.2}
	 * @param versionDir       {@code <catalogue>/<version>/} of the MAG database
	 * @param deleteArchive    delete {@code uhgp-100.tar.gz} after a successful
	 *                         preparation
	 * @return the prepared directory
	 */
	public static File install(String catalogueBaseUrl, File versionDir, boolean deleteArchive, ProgressCallback cb)
			throws IOException {
		ProgressCallback progress = cb == null ? new ProgressCallback.NoOp() : cb;
		File dir = directory(versionDir);
		File archive = new File(dir, PangenomeResource.ARCHIVE);
		if (isInstalled(versionDir)) {
			progress.onMessage("UHGP-100 protein catalogue already prepared: " + dir.getAbsolutePath());
			deleteArchiveIfAsked(archive, deleteArchive, progress);
			return dir;
		}
		if (!dir.isDirectory() && !dir.mkdirs()) {
			throw new IOException("Cannot create " + dir.getAbsolutePath());
		}

		File metadata = new File(versionDir, PangenomeResource.METADATA);
		if (!metadata.isFile()) {
			ResumableDownload.download(catalogueBaseUrl + "/" + PangenomeResource.METADATA, metadata, true, progress);
		}

		if (!archive.isFile()) {
			File partial = new File(dir, PangenomeResource.ARCHIVE + DOWNLOAD_SUFFIX);
			requireSpace(dir, ARCHIVE_BYTES - (partial.isFile() ? partial.length() : 0L)
					+ preparedBytes(ARCHIVE_BYTES));
			progress.onMessage(String.format(Locale.ROOT,
					"Downloading the UHGP-100 protein catalogue (%.1f GB); an interrupted download resumes",
					ARCHIVE_BYTES / 1e9));
			ResumableDownload.download(catalogueBaseUrl + "/" + PangenomeResource.DEFAULT_SUBDIR + "/"
					+ PangenomeResource.ARCHIVE, partial, true, progress);
			Files.move(partial.toPath(), archive.toPath(), StandardCopyOption.REPLACE_EXISTING);
		} else {
			requireSpace(dir, preparedBytes(archive.length()));
		}

		progress.onMessage("Preparing the UHGP-100 protein catalogue (about 45 minutes)...");
		PangenomeResourcePreparer.prepare(archive, metadata, dir, progress::onMessage);
		if (!isInstalled(versionDir)) {
			throw new IOException("Preparation finished but the result in " + dir.getAbsolutePath()
					+ " is incomplete");
		}
		deleteArchiveIfAsked(archive, deleteArchive, progress);
		progress.onMessage("UHGP-100 protein catalogue ready: " + dir.getAbsolutePath());
		return dir;
	}

	private static void deleteArchiveIfAsked(File archive, boolean deleteArchive, ProgressCallback progress) {
		if (deleteArchive && archive.isFile()) {
			if (archive.delete()) {
				progress.onMessage("Deleted " + archive.getName() + " after preparation");
			} else {
				progress.onMessage("WARNING: could not delete " + archive.getAbsolutePath());
			}
		}
	}

	/** The prepared output scales with the archive: 22.9 GB from 28.3 GB for v2.0.2. */
	static long preparedBytes(long archiveBytes) {
		return (long) (archiveBytes * ((double) PREPARED_BYTES / ARCHIVE_BYTES));
	}

	/** {@code bytes} plus a 5% margin. */
	private static void requireSpace(File dir, long bytes) throws IOException {
		long needed = bytes + bytes / 20;
		long usable = dir.getUsableSpace();
		if (usable > 0 && usable < needed) {
			throw new IOException(String.format(Locale.ROOT,
					"Not enough disk space for the UHGP-100 protein catalogue in %s: %.1f GB needed, %.1f GB free",
					dir.getAbsolutePath(), needed / 1e9, usable / 1e9));
		}
	}
}
