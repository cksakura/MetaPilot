package uk.ac.quadram.metapilot.core.variant;

import java.io.BufferedInputStream;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.InterruptedIOException;
import java.io.OutputStream;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.StandardCopyOption;
import java.security.MessageDigest;
import java.security.NoSuchAlgorithmException;
import java.util.BitSet;
import java.util.Locale;
import java.util.function.Consumer;
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;

import org.json.JSONArray;
import org.json.JSONObject;

/**
 * One-time preparation of the UHGP-100 resource from the downloaded
 * {@code uhgp-100.tar.gz} (28.3 GB, MGnify human-gut v2.0.2).
 * <p>
 * The archive holds {@code uhgp-100/uhgp-100.tsv} (cluster membership, 25 GB)
 * and then {@code uhgp-100/uhgp-100.faa} (representative sequences, 41.6 GB),
 * in that order, so one sequential read produces both outputs:
 * <ul>
 * <li>{@code uhgp-100.faa.gz} - the FASTA, re-gzipped, headers untouched. The
 * lookup streams this on every run; it is never decompressed to disk.</li>
 * <li>{@code uhgp-100.rep_species.tsv.gz} - {@code representative \t species,...}
 * for clusters whose members span more than one species. The TSV is not
 * strictly grouped by representative (about 19.6k clusters recur in separate
 * blocks), so rows are unioned per representative before writing.</li>
 * <li>{@code uhgp-100.prepare_manifest.json} - sizes and the SHA-256 of the
 * uncompressed FASTA, so later runs can tell a complete preparation from an
 * interrupted one.</li>
 * </ul>
 * Outputs are written as {@code .part} files and renamed only when everything
 * succeeded. Measured with the equivalent Python script on 2026-10-07: about 43
 * minutes end to end.
 */
public final class PangenomeResourcePreparer {

	private static final int BUFFER = 1 << 20;

	private PangenomeResourcePreparer() {
	}

	/**
	 * @param archive     {@code uhgp-100.tar.gz}
	 * @param metadataTsv {@code genomes-all_metadata.tsv}, or null; without it
	 *                    the species table is not written and the lookup reports
	 *                    the ALL scope only
	 * @param outDir      where the three outputs go
	 */
	public static JSONObject prepare(File archive, File metadataTsv, File outDir, Consumer<String> progress)
			throws IOException {
		outDir.mkdirs();
		File faaPart = new File(outDir, PangenomeResource.FAA + ".part");
		File repPart = new File(outDir, PangenomeResource.REP_SPECIES + ".part");
		boolean done = false;
		try {
			JSONObject manifest = prepare(archive, metadataTsv, outDir, faaPart, repPart,
					progress == null ? s -> {
					} : progress);
			done = true;
			return manifest;
		} finally {
			if (!done) {
				faaPart.delete();
				repPart.delete();
			}
		}
	}

	private static JSONObject prepare(File archive, File metadataTsv, File outDir, File faaPart, File repPart,
			Consumer<String> log) throws IOException {
		SpeciesMap species = metadataTsv != null && metadataTsv.isFile() ? SpeciesMap.metadataOnly(metadataTsv) : null;
		if (species == null) {
			log.accept("Pangenome prepare: no genomes-all_metadata.tsv found; the species table will not be "
					+ "written and only the ALL scope can be reported");
		}

		JSONObject manifest = new JSONObject();
		manifest.put("prepared_by", "MetaPilot " + PangenomeResourcePreparer.class.getSimpleName());
		manifest.put("source_archive", archive.getAbsolutePath());
		manifest.put("source_content_length", archive.length());
		manifest.put("source_last_modified_millis", archive.lastModified());
		JSONArray members = new JSONArray();
		boolean sawFaa = false;
		long start = System.currentTimeMillis();
		try (InputStream raw = new BufferedInputStream(new FileInputStream(archive), BUFFER);
				InputStream gz = new GZIPInputStream(raw, BUFFER)) {
			TarStream tar = new TarStream(gz);
			TarStream.Entry e;
			while ((e = tar.next()) != null) {
				JSONObject m = new JSONObject();
				m.put("name", e.name);
				m.put("size", e.size);
				m.put("type", e.isFile() ? "file" : String.valueOf(e.type));
				members.put(m);
				if (!e.isFile()) {
					continue;
				}
				String lower = e.name.toLowerCase(Locale.ROOT);
				if (lower.endsWith(".tsv") && species != null) {
					log.accept("Pangenome prepare: reading cluster membership " + e.name + " (" + e.size + " bytes)");
					long t = System.currentTimeMillis();
					readMembership(tar.content(), species, manifest, log);
					manifest.put("tsv_seconds", (System.currentTimeMillis() - t) / 1000);
				} else if (lower.endsWith(".faa")) {
					log.accept("Pangenome prepare: re-compressing " + e.name + " (" + e.size + " bytes)");
					long t = System.currentTimeMillis();
					copyFasta(tar.content(), faaPart, manifest, log);
					manifest.put("faa_seconds", (System.currentTimeMillis() - t) / 1000);
					sawFaa = true;
				}
			}
		}
		manifest.put("archive_members", members);
		if (!sawFaa) {
			Files.deleteIfExists(faaPart.toPath());
			throw new IOException("No .faa member found in " + archive);
		}
		if (species != null) {
			int before = species.repRowCount();
			species.sortAndMerge();
			try (OutputStream out = gzip(repPart)) {
				species.writeRepSpecies(out);
			}
			manifest.put("rep_species_rows_before_merge", before);
			manifest.put("rep_species_rows", species.repRowCount());
		}

		File faa = new File(outDir, PangenomeResource.FAA);
		File rep = new File(outDir, PangenomeResource.REP_SPECIES);
		Files.move(faaPart.toPath(), faa.toPath(), StandardCopyOption.REPLACE_EXISTING);
		manifest.put("out_faa_gz_bytes", faa.length());
		if (species != null) {
			Files.move(repPart.toPath(), rep.toPath(), StandardCopyOption.REPLACE_EXISTING);
			manifest.put("out_rep_species_gz_bytes", rep.length());
		}
		Files.write(new File(outDir, PangenomeResource.MANIFEST).toPath(),
				manifest.toString(2).getBytes(StandardCharsets.UTF_8));
		log.accept(String.format(Locale.ROOT, "Pangenome prepare: done in %.0f min; %s written to %s",
				(System.currentTimeMillis() - start) / 60000.0, PangenomeResource.FAA, outDir));
		return manifest;
	}

	/** Usage: {@code <uhgp-100.tar.gz> <genomes-all_metadata.tsv|-> <out_dir>}. */
	public static void main(String[] args) throws IOException {
		if (args.length != 3) {
			System.err.println("Usage: PangenomeResourcePreparer <uhgp-100.tar.gz> <genomes-all_metadata.tsv|-> "
					+ "<out_dir>");
			System.exit(2);
		}
		File metadata = args[1].equals("-") ? null : new File(args[1]);
		prepare(new File(args[0]), metadata, new File(args[2]), System.out::println);
	}

	private static void readMembership(InputStream content, SpeciesMap species, JSONObject manifest,
			Consumer<String> log) throws IOException {
		BufferedReader br = new BufferedReader(new InputStreamReader(content, StandardCharsets.US_ASCII), BUFFER);
		String currentRep = null;
		BitSet clusterSpecies = new BitSet();
		long members = 0;
		long blocks = 0;
		long unknown = 0;
		long unencodable = 0;
		JSONArray unknownExamples = new JSONArray();
		String line;
		while ((line = br.readLine()) != null) {
			int tab = line.indexOf('\t');
			if (tab <= 0) {
				continue;
			}
			members++;
			if (currentRep == null || !line.regionMatches(0, currentRep, 0, tab) || currentRep.length() != tab) {
				if (currentRep != null && !addRow(currentRep, clusterSpecies, species)) {
					unencodable++;
				}
				currentRep = line.substring(0, tab);
				clusterSpecies.clear();
				blocks++;
				if (!addGenomeSpecies(currentRep, clusterSpecies, species)) {
					unknown++;
					if (unknownExamples.length() < 10) {
						unknownExamples.put(currentRep);
					}
				}
			}
			String member = line.substring(tab + 1).trim();
			if (!addGenomeSpecies(member, clusterSpecies, species)) {
				unknown++;
				if (unknownExamples.length() < 10) {
					unknownExamples.put(member);
				}
			}
			if ((members & 0xFFFFF) == 0) {
				checkInterrupted();
			}
			if (members % 100_000_000L == 0) {
				log.accept(String.format(Locale.ROOT, "Pangenome prepare: %,d cluster members read", members));
			}
		}
		if (currentRep != null && !addRow(currentRep, clusterSpecies, species)) {
			unencodable++;
		}
		manifest.put("tsv_members", members);
		manifest.put("tsv_rep_transitions", blocks);
		manifest.put("tsv_unknown_genome_ids", unknown);
		manifest.put("tsv_unknown_genome_examples", unknownExamples);
		manifest.put("tsv_unencodable_multi_species_reps", unencodable);
	}

	private static boolean addGenomeSpecies(String proteinId, BitSet into, SpeciesMap species) {
		int underscore = proteinId.lastIndexOf('_');
		int s = species.speciesOfGenome(underscore > 0 ? proteinId.substring(0, underscore) : proteinId);
		if (s < 0) {
			return false;
		}
		into.set(s);
		return true;
	}

	/** @return false only when a multi-species row could not be stored */
	private static boolean addRow(String rep, BitSet clusterSpecies, SpeciesMap species) {
		if (clusterSpecies.cardinality() <= 1) {
			return true;
		}
		long key = SpeciesMap.proteinKey(rep);
		if (key < 0) {
			return false;
		}
		species.addRepRow(key, clusterSpecies);
		return true;
	}

	private static void copyFasta(InputStream content, File out, JSONObject manifest, Consumer<String> log)
			throws IOException {
		MessageDigest sha;
		try {
			sha = MessageDigest.getInstance("SHA-256");
		} catch (NoSuchAlgorithmException e) {
			throw new IOException(e);
		}
		long records = 0;
		long bytes = 0;
		long nextReport = 5L << 30;
		boolean lineStart = true;
		byte[] buf = new byte[BUFFER];
		try (OutputStream gz = gzip(out)) {
			int n;
			while ((n = content.read(buf)) > 0) {
				checkInterrupted();
				sha.update(buf, 0, n);
				gz.write(buf, 0, n);
				for (int i = 0; i < n; i++) {
					byte b = buf[i];
					if (lineStart && b == '>') {
						records++;
					}
					lineStart = b == '\n';
				}
				bytes += n;
				if (bytes >= nextReport) {
					nextReport += 5L << 30;
					log.accept(String.format(Locale.ROOT, "Pangenome prepare: %,d MB of FASTA, %,d records",
							bytes >> 20, records));
				}
			}
		}
		StringBuilder hex = new StringBuilder();
		for (byte b : sha.digest()) {
			hex.append(String.format(Locale.ROOT, "%02x", b));
		}
		manifest.put("faa_records", records);
		manifest.put("faa_uncompressed_bytes", bytes);
		manifest.put("faa_uncompressed_sha256", hex.toString());
	}

	private static void checkInterrupted() throws InterruptedIOException {
		if (Thread.currentThread().isInterrupted()) {
			throw new InterruptedIOException("Pangenome prepare cancelled");
		}
	}

	private static OutputStream gzip(File f) throws IOException {
		return new GZIPOutputStream(new FileOutputStream(f), BUFFER) {
			{
				def.setLevel(4);
			}
		};
	}
}
