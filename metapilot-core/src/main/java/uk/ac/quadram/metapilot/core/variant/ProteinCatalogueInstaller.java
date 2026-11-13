package uk.ac.quadram.metapilot.core.variant;

import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.net.HttpURLConnection;
import java.net.URI;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.StandardCopyOption;
import java.util.ArrayList;
import java.util.List;
import java.util.Locale;
import java.util.Objects;
import java.util.function.Predicate;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

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
 * Only the human-gut catalogue is supported. The installer reads the
 * {@code protein_catalogue/} directory listing and downloads whichever
 * {@code *-100.tar.gz} is there ({@code uhgp-100.tar.gz} on v2.0.2,
 * {@code protein_catalogue-100.tar.gz} on v3.0, or a later rename of the same
 * 100% archive). It is stored locally as {@code uhgp-100.tar.gz}. The 50/90/95
 * cut-offs and the unclustered {@code proteins.faa.gz} are ignored. Other
 * MGnify catalogues publish protein catalogues under other names.
 */
public final class ProteinCatalogueInstaller {

	public static final String SUPPORTED_CATALOGUE = "human-gut";

	/**
	 * Names probed only when the directory listing itself cannot be read.
	 * A readable listing that has no 100% archive is an error, not a reason to
	 * guess these.
	 */
	static final String[] REMOTE_ARCHIVES = { PangenomeResource.ARCHIVE, "protein_catalogue-100.tar.gz" };

	private static final Pattern HREF = Pattern.compile("(?i)href\\s*=\\s*\"([^\"]+)\"");

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
			File already = existingRemoteArchive(dir);
			if (already != null) {
				progress.onMessage("Using already-downloaded " + already.getName());
				Files.move(already.toPath(), archive.toPath(), StandardCopyOption.REPLACE_EXISTING);
			}
		}
		if (!archive.isFile()) {
			String url = archiveUrl(catalogueBaseUrl);
			long remoteBytes = ResumableDownload.contentLength(url);
			if (remoteBytes <= 0) {
				remoteBytes = ARCHIVE_BYTES;
			}
			File partial = new File(dir, PangenomeResource.ARCHIVE + DOWNLOAD_SUFFIX);
			requireSpace(dir, remoteBytes - (partial.isFile() ? partial.length() : 0L)
					+ preparedBytes(remoteBytes));
			progress.onMessage(String.format(Locale.ROOT,
					"Downloading the UHGP-100 protein catalogue (%.1f GB from %s); an interrupted download resumes",
					remoteBytes / 1e9, remoteName(url)));
			ResumableDownload.download(url, partial, true, progress);
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

	/**
	 * HTTP URL of the 100% protein catalogue archive for this MAG version.
	 * The name comes from the {@code protein_catalogue/} listing.
	 */
	static String archiveUrl(String catalogueBaseUrl) throws IOException {
		String base = catalogueBase(catalogueBaseUrl);
		String directory = base + "/" + PangenomeResource.DEFAULT_SUBDIR + "/";
		String html;
		try {
			html = readDirectoryListing(directory);
		} catch (IOException listingFailed) {
			try {
				return archiveUrl(catalogueBaseUrl, ResumableDownload::remoteExists);
			} catch (IOException probeFailed) {
				listingFailed.addSuppressed(probeFailed);
				throw listingFailed;
			}
		}
		List<String> names = hrefNames(html);
		String chosen = PangenomeResource.selectFullCatalogueArchive(names);
		if (chosen == null) {
			throw new IOException("No 100% protein catalogue archive in " + directory + " (listing: " + names + ")");
		}
		return directory + chosen;
	}

	/** Builds the archive URL from an already-read directory listing. */
	static String archiveUrl(String catalogueBaseUrl, List<String> directoryNames) throws IOException {
		String chosen = PangenomeResource.selectFullCatalogueArchive(directoryNames);
		if (chosen == null) {
			throw new IOException("No 100% protein catalogue archive in listing " + directoryNames);
		}
		return catalogueBase(catalogueBaseUrl) + "/" + PangenomeResource.DEFAULT_SUBDIR + "/" + chosen;
	}

	/**
	 * Probes the historical archive names. Used when the directory listing
	 * cannot be read.
	 */
	static String archiveUrl(String catalogueBaseUrl, Predicate<String> remoteExists) throws IOException {
		Objects.requireNonNull(remoteExists, "remoteExists");
		String base = catalogueBase(catalogueBaseUrl);
		StringBuilder tried = new StringBuilder();
		for (String name : REMOTE_ARCHIVES) {
			String url = base + "/" + PangenomeResource.DEFAULT_SUBDIR + "/" + name;
			if (tried.length() > 0) {
				tried.append(", ");
			}
			tried.append(url);
			if (remoteExists.test(url)) {
				return url;
			}
		}
		throw new IOException("No 100% protein catalogue archive found (tried " + tried + ")");
	}

	private static String catalogueBase(String catalogueBaseUrl) {
		Objects.requireNonNull(catalogueBaseUrl, "catalogueBaseUrl");
		return catalogueBaseUrl.endsWith("/")
				? catalogueBaseUrl.substring(0, catalogueBaseUrl.length() - 1)
				: catalogueBaseUrl;
	}

	/** File names linked from an Apache/nginx directory index. Query links and parent paths are dropped. */
	static List<String> hrefNames(String html) {
		List<String> names = new ArrayList<>();
		if (html == null || html.isBlank()) {
			return names;
		}
		Matcher matcher = HREF.matcher(html);
		while (matcher.find()) {
			String href = matcher.group(1).trim();
			if (href.isEmpty() || href.startsWith("?") || href.startsWith("/") || href.contains("..")
					|| href.endsWith("/")) {
				continue;
			}
			int slash = href.lastIndexOf('/');
			String name = slash >= 0 ? href.substring(slash + 1) : href;
			int query = name.indexOf('?');
			if (query >= 0) {
				name = name.substring(0, query);
			}
			if (!name.isBlank() && !names.contains(name)) {
				names.add(name);
			}
		}
		return names;
	}

	private static String readDirectoryListing(String url) throws IOException {
		String current = url;
		HttpURLConnection http = null;
		for (int redirect = 0; redirect < 5; redirect++) {
			http = (HttpURLConnection) URI.create(current).toURL().openConnection();
			http.setInstanceFollowRedirects(false);
			http.setConnectTimeout(30_000);
			http.setReadTimeout(30_000);
			http.setRequestProperty("User-Agent", "MetaPilot");
			int code = http.getResponseCode();
			if (code >= 300 && code < 400) {
				String location = http.getHeaderField("Location");
				http.disconnect();
				http = null;
				if (location == null || location.isBlank()) {
					throw new IOException("Redirect without a location while listing " + url);
				}
				current = URI.create(current).resolve(location).toString();
				continue;
			}
			if (code < 200 || code >= 300) {
				http.disconnect();
				throw new IOException("HTTP " + code + " listing " + url);
			}
			try (InputStream in = http.getInputStream()) {
				byte[] body = in.readNBytes(1_000_000);
				return new String(body, StandardCharsets.ISO_8859_1);
			} finally {
				http.disconnect();
			}
		}
		if (http != null) {
			http.disconnect();
		}
		throw new IOException("Too many redirects while listing " + url);
	}

	private static String remoteName(String url) {
		int slash = url.lastIndexOf('/');
		return slash >= 0 ? url.substring(slash + 1) : url;
	}

	/** A complete 100% archive already on disk under a name other than {@link PangenomeResource#ARCHIVE}. */
	private static File existingRemoteArchive(File dir) {
		File found = PangenomeResource.fullCatalogueArchive(dir);
		if (found == null || PangenomeResource.ARCHIVE.equalsIgnoreCase(found.getName())) {
			return null;
		}
		return found;
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
