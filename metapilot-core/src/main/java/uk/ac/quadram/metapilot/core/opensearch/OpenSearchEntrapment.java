package uk.ac.quadram.metapilot.core.opensearch;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.StandardCopyOption;
import java.security.MessageDigest;
import java.security.NoSuchAlgorithmException;
import java.util.HashMap;
import java.util.HashSet;
import java.util.HexFormat;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.zip.GZIPInputStream;

import uk.ac.quadram.metapilot.core.dia.EntrapmentShuffle;
import uk.ac.quadram.metapilot.core.prodb.TargetDecoyFasta;

/**
 * Evaluation-only, round-two entrapment. Controls remain targets during identification.
 * CLI ({@code mag-opensearch-entrapment}) and the MAG open-search GUI switch both use this path.
 */
public final class OpenSearchEntrapment {
    public static final String RESULT_FOLDER = "opensearch_entrapment_MAG_search";
    public static final long DEFAULT_SEED = 1L;

    private OpenSearchEntrapment() { }

    public static final class Config {
        public final long seed;
        public final int idDigit;

        public Config(long seed, int idDigit) {
            if (idDigit < 1 || idDigit > 9) throw new IllegalArgumentException("Entrapment ID digit must be 1-9");
            this.seed = seed;
            this.idDigit = idDigit;
        }
    }

    /**
     * Keeps the real-only FASTA intact. Rebuilds the searched FASTA from real + shuffled
     * targets, then generates rev_ decoys of both. The last-written manifest commits
     * the source digest, seed, ID namespace and generated digests for safe resume.
     */
    public static void prepare(File realTargets, File searchedFasta, Config config) throws IOException {
        prepare(realTargets, searchedFasta, config, TargetDecoyFasta.DecoyMethod.REVERSE);
    }

    /**
     * As {@link #prepare(File, File, Config)}, with the decoys of both populations made by
     * {@code decoyMethod}. A non-default method is part of the manifest identity, so switching it
     * rebuilds the searched FASTA; the default leaves existing manifests valid.
     */
    public static void prepare(File realTargets, File searchedFasta, Config config,
            TargetDecoyFasta.DecoyMethod decoyMethod) throws IOException {
        Path source = realTargets.toPath();
        File augmented = new File(realTargets.getParentFile(), realTargets.getName() + ".entrapment.fasta");
        Path manifest = Path.of(searchedFasta.getPath() + ".entrapment.tsv");
        String identity = "format\t1\nseed\t" + config.seed + "\nid_digit\t" + config.idDigit
                + "\nreal_targets_sha256\t" + digest(source) + "\n"
                + (decoyMethod == TargetDecoyFasta.DecoyMethod.REVERSE ? "" : "decoy_method\t" + decoyMethod.label + "\n");
        if (Files.isRegularFile(manifest) && augmented.isFile() && searchedFasta.isFile()) {
            String old = Files.readString(manifest);
            if (old.startsWith(identity)
                    && old.contains("augmented_sha256\t" + digest(augmented.toPath()) + "\n")
                    && old.contains("searched_sha256\t" + digest(searchedFasta.toPath()) + "\n")
                    && TargetDecoyFasta.isValid(augmented, searchedFasta, decoyMethod)) return;
        }
        // Reject decoys, existing controls and duplicate IDs before constructing the paired space.
        Set<String> ids = new HashSet<>();
        try (BufferedReader reader = Files.newBufferedReader(source)) {
            String line;
            while ((line = reader.readLine()) != null) {
                if (!line.startsWith(">")) continue;
                String id = line.substring(1).trim().split("\\s+", 2)[0];
                if (!ids.add(id)) throw new IOException("Duplicate real protein ID: " + id);
                try { EntrapmentShuffle.entrapmentId(id, config.idDigit); }
                catch (IllegalArgumentException ex) {
                    throw new IOException("Entrapment requires real MGYG0 protein IDs; found " + id, ex);
                }
            }
        }
        if (ids.isEmpty()) throw new IOException("No real proteins to pair in " + realTargets);
        Path temp = Files.createTempFile(realTargets.toPath().getParent(), "entrapment-", ".fasta");
        try {
            int shuffled;
            try (PrintWriter writer = new PrintWriter(Files.newBufferedWriter(temp, StandardCharsets.UTF_8))) {
                try (BufferedReader reader = Files.newBufferedReader(source)) {
                    String line;
                    while ((line = reader.readLine()) != null) writer.println(line);
                }
                shuffled = EntrapmentShuffle.appendShuffledGenome(realTargets, writer, config.idDigit, config.seed);
                if (writer.checkError()) throw new IOException("Could not write " + temp);
            }
            if (shuffled != ids.size()) throw new IOException("Empty or missing sequences in paired FASTA: "
                    + shuffled + " shuffled proteins for " + ids.size() + " real headers");
            // A failed rebuild cannot leave the old manifest certifying a partially replaced database.
            Files.deleteIfExists(manifest);
            Files.move(temp, augmented.toPath(), StandardCopyOption.REPLACE_EXISTING);
            // Explicit invalidation also covers seed changes on coarse-resolution filesystems.
            Files.deleteIfExists(searchedFasta.toPath());
            TargetDecoyFasta.ensure(augmented, searchedFasta, decoyMethod);
            String record = identity + "real_proteins\t" + ids.size() + "\nshuffled_proteins\t" + shuffled
                    + "\nreversed_decoys\t" + (2L * ids.size())
                    + "\naugmented_sha256\t" + digest(augmented.toPath())
                    + "\nsearched_sha256\t" + digest(searchedFasta.toPath()) + "\n";
            Files.writeString(manifest, record, StandardCharsets.UTF_8);
        } finally {
            Files.deleteIfExists(temp);
        }
    }

    private static String digest(Path file) throws IOException {
        try {
            MessageDigest digest = MessageDigest.getInstance("SHA-256");
            try (InputStream in = Files.newInputStream(file)) {
                byte[] buf = new byte[1 << 16];
                int n;
                while ((n = in.read(buf)) != -1) digest.update(buf, 0, n);
            }
            return HexFormat.of().formatHex(digest.digest());
        } catch (NoSuchAlgorithmException ex) { throw new IllegalStateException(ex); }
    }

    /**
     * A descriptive accession-based exposure readout, not an FDP estimate. The arm
     * map is search-reported, not exhaustive sequence remapping; shared/unshufflable
     * peptides need independent collision handling before estimating error rates.
     */
    public static void writeSummary(File combined, Config config) throws IOException {
        Map<String, Integer> labels = new HashMap<>();
        try (BufferedReader reader = gzip(new File(combined, "identification_arm_peptides.tsv.gz"))) {
            reader.readLine();
            String line;
            while ((line = reader.readLine()) != null) {
                String[] c = line.split("\t", -1);
                int flags = 0;
                for (String id : c[1].split(";")) {
                    if (id.isBlank() || EntrapmentShuffle.isDecoyId(id)) continue;
                    flags |= EntrapmentShuffle.isEntrapment(id, config.idDigit) ? 2 : 1;
                }
                labels.merge(base(c[0]), flags, (a, b) -> a | b);
            }
        }
        Map<String, long[]> counts = new TreeMap<>();
        Map<String, Set<String>> peptides = new HashMap<>();
        try (BufferedReader reader = gzip(new File(combined, "identification_arms.tsv.gz"))) {
            reader.readLine();
            String line;
            while ((line = reader.readLine()) != null) {
                String[] c = line.split("\t", -1);
                if (!Boolean.parseBoolean(c[7])) continue;
                String label = !Boolean.parseBoolean(c[4]) ? "decoy" : switch (labels.getOrDefault(base(c[3]), 0)) {
                    case 1 -> "real_accessions";
                    case 2 -> "entrapment_accessions";
                    case 3 -> "shared_accessions";
                    default -> "unmapped";
                };
                String key = c[0] + "\t" + label;
                long[] n = counts.computeIfAbsent(key, k -> new long[1]);
                n[0]++;
                peptides.computeIfAbsent(key, k -> new HashSet<>()).add(base(c[3]));
            }
        }
        try (PrintWriter writer = new PrintWriter(new File(combined, "entrapment_summary.tsv"), StandardCharsets.UTF_8)) {
            writer.println("# seed=" + config.seed + "; id_digit=" + config.idDigit + "; 1 shuffled protein per real protein");
            writer.println("# Saved arm acceptance flags; real and entrapment compete as targets, each with reversed decoys.");
            writer.println("# Descriptive accession labels only, NOT an FDP estimate. Search-reported protein mappings");
            writer.println("# may omit real/shuffled peptide collisions. Resolve collisions before measuring identification error.");
            writer.println("# Arms retain their own acceptance rules; counts are not a matched-error comparison.");
            writer.println("arm\tpopulation\taccepted_psms\tdistinct_base_peptides_IL_equivalent");
            for (Map.Entry<String, long[]> e : counts.entrySet()) {
                writer.println(e.getKey() + "\t" + e.getValue()[0] + "\t" + peptides.get(e.getKey()).size());
            }
            if (writer.checkError()) throw new IOException("Could not write entrapment summary in " + combined);
        }
    }

    private static String base(String peptide) { return peptide.toUpperCase(java.util.Locale.ROOT).replace('I', 'L'); }

    private static BufferedReader gzip(File file) throws IOException {
        return new BufferedReader(new InputStreamReader(new GZIPInputStream(Files.newInputStream(file.toPath())),
                StandardCharsets.UTF_8));
    }
}
