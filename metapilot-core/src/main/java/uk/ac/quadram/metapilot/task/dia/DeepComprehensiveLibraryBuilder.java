package uk.ac.quadram.metapilot.task.dia;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import uk.ac.quadram.metapilot.task.mag.MagDbItem;

/**
 * Builds a FASTA file containing only the requested protein IDs from the MAG DB.
 * Used for iterative Deep Comprehensive rounds: library for round N contains all
 * peptides identified in the previous round (via their proteins).
 */
public final class DeepComprehensiveLibraryBuilder {

	private static final Logger LOG = LogManager.getLogger(DeepComprehensiveLibraryBuilder.class);

	private DeepComprehensiveLibraryBuilder() {}

	/** Extract genome ID from protein ID (prefix before first '_'). */
	public static String genomeFromProteinId(String proteinId) {
		if (proteinId == null) return "";
		int u = proteinId.indexOf('_');
		return u > 0 ? proteinId.substring(0, u) : proteinId;
	}

	/**
	 * Keep only protein IDs whose genome exists in the MAG DB (genome.faa or genome.fasta).
	 * Use this before writeFastaFromProteinIds when the set may contain IDs from non-MAG libraries
	 * (e.g. user library like "pro17511"); those are dropped so we don't log "Genome file not found"
	 * for every one and get 0 entries written.
	 *
	 * @param proteinIds set of protein IDs (may mix MAG-style genome_proteinId and other names)
	 * @param magDb      MAG DB (getOriginalDbFolder())
	 * @return new set containing only IDs whose genome file exists in the MAG DB
	 */
	public static Set<String> filterProteinIdsToMagGenomes(Set<String> proteinIds, MagDbItem magDb) {
		Set<String> out = new HashSet<>();
		if (proteinIds == null || proteinIds.isEmpty() || magDb == null) return out;
		File dbFolder = magDb.getOriginalDbFolder();
		if (dbFolder == null || !dbFolder.isDirectory()) return out;
		for (String pid : proteinIds) {
			String genome = genomeFromProteinId(pid);
			if (genome.isEmpty()) continue;
			File genomeFaa = new File(dbFolder, genome + ".faa");
			if (!genomeFaa.isFile()) {
				genomeFaa = new File(dbFolder, genome + ".fasta");
			}
			if (genomeFaa.isFile()) {
				out.add(pid);
			}
		}
		return out;
	}

	/**
	 * Write a FASTA file with only the entries for the given protein IDs.
	 * Reads from MAG DB original folder (genome.faa or genome.fasta per genome).
	 *
	 * @param proteinIds set of full protein IDs (e.g. MGYG000000001_00001)
	 * @param magDb      MAG DB (getOriginalDbFolder())
	 * @param outputFasta output FASTA path
	 * @return number of protein entries written
	 */
	public static int writeFastaFromProteinIds(Set<String> proteinIds, MagDbItem magDb, File outputFasta) throws IOException {
		if (proteinIds == null || proteinIds.isEmpty()) {
			LOG.warn("No protein IDs provided");
			return 0;
		}
		if (magDb == null) {
			throw new IOException("MagDbItem is required");
		}
		File dbFolder = magDb.getOriginalDbFolder();
		if (dbFolder == null || !dbFolder.isDirectory()) {
			throw new IOException("MAG DB original folder not available");
		}
		// Group protein IDs by genome
		Map<String, Set<String>> genomeToProteins = new HashMap<>();
		for (String pid : proteinIds) {
			String genome = genomeFromProteinId(pid);
			if (genome.isEmpty()) continue;
			genomeToProteins.computeIfAbsent(genome, k -> new HashSet<>()).add(pid);
		}
		File parent = outputFasta.getParentFile();
		if (parent != null && !parent.exists()) {
			parent.mkdirs();
		}
		int written = 0;
		try (PrintWriter w = new PrintWriter(outputFasta, "UTF-8")) {
			for (Map.Entry<String, Set<String>> e : genomeToProteins.entrySet()) {
				String genome = e.getKey();
				Set<String> want = e.getValue();
				File genomeFaa = new File(dbFolder, genome + ".faa");
				if (!genomeFaa.exists()) {
					genomeFaa = new File(dbFolder, genome + ".fasta");
				}
				if (!genomeFaa.isFile()) {
					LOG.warn("Genome file not found: {}", genome);
					continue;
				}
				try (BufferedReader r = new BufferedReader(new FileReader(genomeFaa))) {
					String line;
					boolean inEntry = false;
					String currentId = null;
					while ((line = r.readLine()) != null) {
						if (line.startsWith(">")) {
							if (inEntry && currentId != null && want.contains(currentId)) {
								written++;
							}
							String header = line.substring(1).trim();
							int firstSpace = header.indexOf(' ');
							currentId = firstSpace > 0 ? header.substring(0, firstSpace).trim() : header.trim();
							inEntry = want.contains(currentId);
							if (inEntry) {
								w.println(line);
							}
						} else if (inEntry) {
							w.println(line);
						}
					}
					if (inEntry && currentId != null && want.contains(currentId)) {
						written++;
					}
				}
			}
		}
		LOG.info("Wrote {} protein entries to {}", written, outputFasta.getAbsolutePath());
		return written;
	}
}
