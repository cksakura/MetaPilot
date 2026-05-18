import java.io.IOException;
import java.nio.file.DirectoryStream;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.StandardCopyOption;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;

public class FoldNameUpdate {

	private static final String OLD_MARKER_FOLDER = "hap";
	private static final String NEW_MARKER_FOLDER = "marker";
	private static final String NEW_SAMPLE_DB_FOLDER = "sample_db";
	private static final String[] OLD_SAMPLE_DB_FOLDERS = { "sample-db", "samaple-db", "sampledb" };
	private static final String[] SEARCH_RESULT_FILES = {
			"results.sage.tsv",
			"pFind.cfg",
			"pFind.protein",
			"pFind.spectra",
			"pFind.summary",
			"pFind-Filtered.spectra"
	};
	private static final String[] OLD_MARKER_PREFIX_FILES = {
			"marker_genome_evidence.tsv",
			"marker_genome_ranked.tsv",
			"marker_genome_selected.tsv"
	};

	public static void main(String[] args) throws IOException {
		if (args.length == 0) {
			System.out.println("Usage: java FoldNameUpdate <result-folder> [--dry-run]");
			System.out.println("Example: java FoldNameUpdate Z:\\PRIDE_dataset\\PXD...\\MetaPilot --dry-run");
			return;
		}

		Path root = Path.of(args[0]);
		boolean dryRun = args.length > 1 && "--dry-run".equalsIgnoreCase(args[1]);
		UpdateSummary summary = updateFolderNames(root, dryRun);
		System.out.println(summary);
	}

	/**
	 * Updates legacy result folders so current post-analysis code can read them.
	 *
	 * <ul>
	 * <li>hap -> marker</li>
	 * <li>create sample_db and copy per-raw search results into it if missing</li>
	 * <li>copy sample-db / samaple-db / sampledb contents into sample_db</li>
	 * <li>marker_genome_*.tsv -> genome_*.tsv aliases when the new names are absent</li>
	 * </ul>
	 */
	public static UpdateSummary updateFolderNames(Path resultRoot, boolean dryRun) throws IOException {
		if (resultRoot == null) {
			throw new IllegalArgumentException("resultRoot cannot be null");
		}
		if (!Files.isDirectory(resultRoot)) {
			throw new IOException("Result folder is not a directory: " + resultRoot.toAbsolutePath());
		}

		UpdateSummary summary = new UpdateSummary();
		List<Path> directories = listDirectoriesDeepestFirst(resultRoot);
		for (Path directory : directories) {
			Path current = directory;
			if (!Files.isDirectory(current)) {
				continue;
			}

			String name = current.getFileName().toString();
			if (OLD_MARKER_FOLDER.equals(name)) {
				current = moveOrMergeDirectory(current, current.resolveSibling(NEW_MARKER_FOLDER), dryRun, summary);
			} else if (isOldSampleDbFolder(name)) {
				copyDirectoryContents(current, current.resolveSibling(NEW_SAMPLE_DB_FOLDER), dryRun, summary);
			}

			if (Files.isDirectory(current) && NEW_MARKER_FOLDER.equals(current.getFileName().toString())) {
				normalizeMarkerSelectionFileNames(current, dryRun, summary);
			}
			if (Files.isDirectory(current)) {
				copyRawSearchResultsToSampleDb(current, dryRun, summary);
			}
		}

		return summary;
	}

	private static List<Path> listDirectoriesDeepestFirst(Path root) throws IOException {
		List<Path> directories = new ArrayList<>();
		try (var stream = Files.walk(root)) {
			stream.filter(Files::isDirectory).forEach(directories::add);
		}
		directories.sort(Comparator.comparingInt(Path::getNameCount).reversed());
		return directories;
	}

	private static boolean isOldSampleDbFolder(String name) {
		for (String oldName : OLD_SAMPLE_DB_FOLDERS) {
			if (oldName.equalsIgnoreCase(name)) {
				return true;
			}
		}
		return false;
	}

	private static Path moveOrMergeDirectory(Path source, Path target, boolean dryRun, UpdateSummary summary)
			throws IOException {
		if (source.equals(target)) {
			return source;
		}

		if (!Files.exists(target)) {
			summary.recordRename(source, target);
			if (!dryRun) {
				Files.move(source, target);
			}
			return target;
		}

		if (!Files.isDirectory(target)) {
			summary.recordSkip(source, "target exists and is not a directory: " + target);
			return source;
		}

		summary.recordMerge(source, target);
		if (!dryRun) {
			moveMissingDirectoryContents(source, target, summary);
		}
		return target;
	}

	private static void moveMissingDirectoryContents(Path source, Path target, UpdateSummary summary) throws IOException {
		try (DirectoryStream<Path> entries = Files.newDirectoryStream(source)) {
			for (Path sourceChild : entries) {
				Path targetChild = target.resolve(sourceChild.getFileName().toString());
				if (!Files.exists(targetChild)) {
					Files.move(sourceChild, targetChild);
					continue;
				}
				if (Files.isDirectory(sourceChild) && Files.isDirectory(targetChild)) {
					moveMissingDirectoryContents(sourceChild, targetChild, summary);
				} else {
					summary.recordSkip(sourceChild, "target already exists: " + targetChild);
				}
			}
		}
	}

	private static void copyDirectoryContents(Path source, Path target, boolean dryRun, UpdateSummary summary)
			throws IOException {
		if (!Files.isDirectory(source)) {
			return;
		}
		if (!Files.exists(target)) {
			summary.recordCreate(target);
			if (!dryRun) {
				Files.createDirectories(target);
			}
		}

		try (DirectoryStream<Path> entries = Files.newDirectoryStream(source)) {
			for (Path sourceChild : entries) {
				Path targetChild = target.resolve(sourceChild.getFileName().toString());
				if (Files.isDirectory(sourceChild)) {
					copyDirectoryContents(sourceChild, targetChild, dryRun, summary);
				} else if (!Files.exists(targetChild)) {
					summary.recordCopy(sourceChild, targetChild);
					if (!dryRun) {
						Files.copy(sourceChild, targetChild, StandardCopyOption.COPY_ATTRIBUTES);
					}
				} else {
					summary.recordSkip(sourceChild, "target already exists: " + targetChild);
				}
			}
		}
	}

	private static void copyRawSearchResultsToSampleDb(Path rawFolder, boolean dryRun, UpdateSummary summary)
			throws IOException {
		if (NEW_SAMPLE_DB_FOLDER.equals(rawFolder.getFileName().toString())
				|| NEW_MARKER_FOLDER.equals(rawFolder.getFileName().toString())
				|| OLD_MARKER_FOLDER.equals(rawFolder.getFileName().toString())
				|| isOldSampleDbFolder(rawFolder.getFileName().toString())) {
			return;
		}

		Path sampleDb = rawFolder.resolve(NEW_SAMPLE_DB_FOLDER);
		for (String resultFileName : SEARCH_RESULT_FILES) {
			Path source = rawFolder.resolve(resultFileName);
			if (Files.isRegularFile(source)) {
				copyFileToSampleDb(source, sampleDb.resolve(resultFileName), dryRun, summary);
			}
		}

		Path marker = rawFolder.resolve(NEW_MARKER_FOLDER);
		if (Files.isDirectory(marker)) {
			for (String resultFileName : SEARCH_RESULT_FILES) {
				Path source = marker.resolve(resultFileName);
				if (Files.isRegularFile(source)) {
					copyFileToSampleDb(source, sampleDb.resolve(resultFileName), dryRun, summary);
				}
			}
		}
	}

	private static void copyFileToSampleDb(Path source, Path target, boolean dryRun, UpdateSummary summary)
			throws IOException {
		if (Files.exists(target)) {
			summary.recordSkip(source, "target already exists: " + target);
			return;
		}
		Path parent = target.getParent();
		if (!Files.exists(parent)) {
			summary.recordCreate(parent);
			if (!dryRun) {
				Files.createDirectories(parent);
			}
		}

		summary.recordCopy(source, target);
		if (!dryRun) {
			Files.copy(source, target, StandardCopyOption.COPY_ATTRIBUTES);
		}
	}

	private static void normalizeMarkerSelectionFileNames(Path markerFolder, boolean dryRun, UpdateSummary summary)
			throws IOException {
		for (String oldFileName : OLD_MARKER_PREFIX_FILES) {
			Path oldFile = markerFolder.resolve(oldFileName);
			if (!Files.isRegularFile(oldFile)) {
				continue;
			}

			String newFileName = oldFileName.substring("marker_".length());
			Path newFile = markerFolder.resolve(newFileName);
			if (Files.exists(newFile)) {
				summary.recordSkip(oldFile, "target already exists: " + newFile);
				continue;
			}

			summary.recordCopy(oldFile, newFile);
			if (!dryRun) {
				Files.copy(oldFile, newFile, StandardCopyOption.COPY_ATTRIBUTES);
			}
		}
	}

	public static class UpdateSummary {
		private int renamed;
		private int merged;
		private int created;
		private int copied;
		private int skipped;

		private void recordRename(Path source, Path target) {
			renamed++;
			System.out.println("Rename: " + source + " -> " + target);
		}

		private void recordMerge(Path source, Path target) {
			merged++;
			System.out.println("Merge: " + source + " -> " + target);
		}

		private void recordCreate(Path target) {
			created++;
			System.out.println("Create: " + target);
		}

		private void recordCopy(Path source, Path target) {
			copied++;
			System.out.println("Copy: " + source + " -> " + target);
		}

		private void recordSkip(Path source, String reason) {
			skipped++;
			System.out.println("Skip: " + source + " (" + reason + ")");
		}

		@Override
		public String toString() {
			return "UpdateSummary{renamed=" + renamed
					+ ", merged=" + merged
					+ ", created=" + created
					+ ", copied=" + copied
					+ ", skipped=" + skipped + "}";
		}
	}
}
