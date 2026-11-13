package uk.ac.quadram.metapilot.core.runner;

import uk.ac.quadram.metapilot.core.prodb.TargetDecoyFasta;
import uk.ac.quadram.metapilot.task.par.MetaBaseParameter;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.StandardCopyOption;

/**
 * The host proteome ({@code host_database} with {@code append_host_db}) in the open-search
 * workflow's stage-2 databases.
 * <p>
 * The host joins the sample-specific database that stage 2 and its closed stage search, as the
 * closed FragPipe workflow adds it to its final search: after genome selection, with a decoy per
 * host protein made the same way as the genomes' decoys. It is not added to the marker search,
 * which only selects genomes. Before 2026-11-03 the open-search runner ignored both settings, so
 * runs that set them searched the genomes alone.
 * <p>
 * The genome-only targets stay in {@code <sample>_refined_targets.fasta}; the host is appended in
 * a second file, {@code <sample>_refined_targets_with_host.fasta}, so switching the host on or off
 * rebuilds the searched database without repeating genome selection.
 */
public final class OpenSearchHostDatabase {

	/** Beside the combined targets: the host file the combined targets were built from. */
	static final String STAMP_SUFFIX = ".host";

	private OpenSearchHostDatabase() {
	}

	/**
	 * @return the host FASTA, or null when {@code append_host_db} is off
	 * @throws IOException when the host is switched on but the file is not given or not found:
	 *                     a run that asked for the host must not search without it
	 */
	public static File resolve(MetaBaseParameter params) throws IOException {
		if (params == null || !params.isAppendHostDb()) {
			return null;
		}
		String path = params.getHostDB();
		if (path == null || path.isBlank()) {
			throw new IOException("append_host_db is set but host_database is empty");
		}
		File host = resolveUserPath(path);
		if (!host.isFile() || host.length() == 0) {
			throw new IOException("append_host_db is set but the host database was not found: "
					+ host.getAbsolutePath());
		}
		return host;
	}

	/**
	 * The targets stage 2 searches: {@code genomeTargets} alone without a host, otherwise the
	 * genome targets followed by the host proteins. The combined file is rewritten only when it is
	 * missing, older than either source, or was built from another host file, so an unchanged
	 * setting leaves it - and the database and searches derived from it - as they are.
	 */
	public static File searchTargets(File genomeTargets, File host) throws IOException {
		if (host == null) {
			return genomeTargets;
		}
		String name = genomeTargets.getName();
		String base = name.endsWith(".fasta") ? name.substring(0, name.length() - ".fasta".length()) : name;
		File combined = new File(genomeTargets.getParentFile(), base + "_with_host.fasta");
		File stamp = new File(combined.getPath() + STAMP_SUFFIX);
		String identity = host.getAbsolutePath() + "\t" + host.length() + "\t" + host.lastModified();
		if (combined.isFile() && stamp.isFile() && combined.lastModified() >= genomeTargets.lastModified()
				&& combined.lastModified() >= host.lastModified()
				&& identity.equals(Files.readString(stamp.toPath(), StandardCharsets.UTF_8).trim())) {
			return combined;
		}
		File temp = new File(combined.getParentFile(), combined.getName() + ".tmp." + System.nanoTime());
		try {
			try (PrintWriter writer = new PrintWriter(temp, StandardCharsets.UTF_8)) {
				try (BufferedReader reader = Files.newBufferedReader(genomeTargets.toPath(), StandardCharsets.UTF_8)) {
					String line;
					while ((line = reader.readLine()) != null) {
						writer.println(line);
					}
				}
				int written = appendHost(host, writer);
				if (written == 0) {
					throw new IOException("Host database has no protein records: " + host.getAbsolutePath());
				}
				if (writer.checkError()) {
					throw new IOException("Could not write " + temp);
				}
			}
			Files.move(temp.toPath(), combined.toPath(), StandardCopyOption.REPLACE_EXISTING);
		} finally {
			Files.deleteIfExists(temp.toPath());
		}
		Files.writeString(stamp.toPath(), identity + "\n", StandardCharsets.UTF_8);
		return combined;
	}

	/**
	 * Copies the host's target records. Records already carrying the decoy tag are left out: the
	 * searched database gets exactly one decoy per target from {@link TargetDecoyFasta}.
	 *
	 * @return the number of records written
	 */
	static int appendHost(File host, PrintWriter writer) throws IOException {
		int written = 0;
		boolean keep = false;
		try (BufferedReader reader = Files.newBufferedReader(host.toPath(), StandardCharsets.UTF_8)) {
			String line;
			while ((line = reader.readLine()) != null) {
				line = line.trim();
				if (line.isEmpty()) {
					continue;
				}
				if (line.startsWith(">")) {
					keep = !line.substring(1).startsWith(TargetDecoyFasta.DECOY_TAG);
					if (keep) {
						written++;
					}
				}
				if (keep) {
					writer.println(line);
				}
			}
		}
		return written;
	}

	/**
	 * Whether {@code accession} (a target, or a decoy with the tag removed) is a host protein:
	 * anything that is not a MAG protein of {@code magIdentifier}.
	 */
	public static boolean isHostAccession(String accession, String magIdentifier) {
		return accession != null && magIdentifier != null && !magIdentifier.isEmpty()
				&& !accession.startsWith(magIdentifier);
	}

	private static File resolveUserPath(String path) {
		String trimmed = path.trim();
		if (trimmed.equals("~") || trimmed.startsWith("~/") || trimmed.startsWith("~\\")) {
			String home = System.getProperty("user.home", "");
			if (!home.isBlank()) {
				return trimmed.length() <= 2 ? new File(home) : new File(home, trimmed.substring(2));
			}
		}
		return new File(trimmed);
	}
}
