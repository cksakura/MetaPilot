import java.io.File;
import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Collections;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.StandardCopyOption;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.stream.Stream;

public class MetaPilotDDADIATest {

    /** Equate I and L for peptide comparisons (map I -> L). */
    private static String normalizeIL(String peptide) {
        if (peptide == null || peptide.isEmpty()) return "";
        return peptide.replace('I', 'L');
    }

    /**
     * Re-generate legacy taxonomy_analysis/Taxa.tsv from an existing taxonomy_analysis/Genome.tsv.
     *
     * Legacy behavior (matches the previous exporter):
     * - DIA-style Genome.tsv (has column "S"): per-taxon score = product(genome S), exported as -log10(product).
     * - DDA-style Genome.tsv (has column "Minus_log(p-value)" or "Score"): per-taxon score = sum(genome score),
     *   exported directly.
     *
     * Taxonomy source: GTDB columns in Genome.tsv (GTDB_* ranks).
     *
     * @param genomeTsv Input Genome.tsv (from taxonomy_analysis/Genome.tsv)
     * @param outTaxaTsv Output Taxa.tsv path
     */
    public static void exportLegacyTaxaFromGenomeTsv(File genomeTsv, File outTaxaTsv) throws IOException {
        exportLegacyTaxaFromGenomeTsv(genomeTsv, outTaxaTsv, true);
    }

    /**
     * @param legacyDiaProduct If true, uses legacy DIA aggregation: product(S) then -log10(product).
     *                         If false, uses sum(S) (new-style) for DIA-mode as well.
     */
    public static void exportLegacyTaxaFromGenomeTsv(File genomeTsv, File outTaxaTsv, boolean legacyDiaProduct)
            throws IOException {
        if (genomeTsv == null || !genomeTsv.isFile()) {
            throw new IOException("Genome.tsv not found: " + (genomeTsv == null ? "null" : genomeTsv.getAbsolutePath()));
        }
        File parent = outTaxaTsv.getParentFile();
        if (parent != null && !parent.exists() && !parent.mkdirs()) {
            throw new IOException("Cannot create output folder: " + parent.getAbsolutePath());
        }

        try (BufferedReader reader = new BufferedReader(new FileReader(genomeTsv));
                PrintWriter taxaWriter = new PrintWriter(outTaxaTsv)) {
            String header = reader.readLine();
            if (header == null) {
                throw new IOException("Genome.tsv is empty: " + genomeTsv.getAbsolutePath());
            }
            String[] h = header.split("\t", -1);

            int sIdx = indexOf(h, "S");
            boolean isDia = sIdx >= 0;
            int ddaScoreIdx = -1;
            if (!isDia) {
                ddaScoreIdx = indexOf(h, "Minus_log(p-value)");
                if (ddaScoreIdx < 0) {
                    ddaScoreIdx = indexOf(h, "Score");
                }
                if (ddaScoreIdx < 0) {
                    throw new IOException("Cannot detect genome score column in Genome.tsv header (missing S / Minus_log(p-value) / Score)");
                }
            }

            // GTDB rank columns in order (Domain..Species, etc. as exported)
            ArrayList<Integer> gtdbRankIdx = new ArrayList<>();
            ArrayList<String> rankNames = new ArrayList<>();
            for (int i = 0; i < h.length; i++) {
                String col = h[i].trim();
                if (col.startsWith("GTDB_") && !col.equals("GTDB_Name") && !col.equals("GTDB_Rank")) {
                    gtdbRankIdx.add(i);
                    rankNames.add(col.substring("GTDB_".length()));
                }
            }
            if (gtdbRankIdx.isEmpty()) {
                throw new IOException("No GTDB_* rank columns found in Genome.tsv header: " + genomeTsv.getAbsolutePath());
            }

            // Intensity columns in order
            ArrayList<Integer> intensityIdx = new ArrayList<>();
            ArrayList<String> intensityNames = new ArrayList<>();
            for (int i = 0; i < h.length; i++) {
                String col = h[i];
                if (col != null && col.startsWith("Intensity ")) {
                    intensityIdx.add(i);
                    intensityNames.add(col);
                }
            }

            @SuppressWarnings("unchecked")
            HashMap<String, String[]>[] taxRankMaps = new HashMap[gtdbRankIdx.size()];
            @SuppressWarnings("unchecked")
            HashMap<String, double[]>[] taxIntensityMaps = new HashMap[gtdbRankIdx.size()];
            @SuppressWarnings("unchecked")
            HashMap<String, Double>[] taxScoreMaps = new HashMap[gtdbRankIdx.size()];
            for (int i = 0; i < gtdbRankIdx.size(); i++) {
                taxRankMaps[i] = new HashMap<>();
                taxIntensityMaps[i] = new HashMap<>();
                taxScoreMaps[i] = new HashMap<>();
            }

            String line;
            while ((line = reader.readLine()) != null) {
                if (line.trim().isEmpty()) continue;
                String[] parts = line.split("\t", -1);

                double genomeScore;
                if (isDia) {
                    genomeScore = parseDoubleSafe(parts, sIdx, 0.0d);
                } else {
                    genomeScore = parseDoubleSafe(parts, ddaScoreIdx, 0.0d);
                }

                double[] genomeIntensities = new double[intensityIdx.size()];
                for (int k = 0; k < intensityIdx.size(); k++) {
                    genomeIntensities[k] = parseDoubleSafe(parts, intensityIdx.get(k), 0.0d);
                }

                for (int r = 0; r < gtdbRankIdx.size(); r++) {
                    String taxon = get(parts, gtdbRankIdx.get(r)).trim();
                    if (taxon.isEmpty()) continue;

                    // lineage = all GTDB ranks from 0..r (inclusive)
                    if (!taxRankMaps[r].containsKey(taxon)) {
                        String[] lineage = new String[r + 1];
                        for (int j = 0; j <= r; j++) {
                            lineage[j] = get(parts, gtdbRankIdx.get(j)).trim();
                        }
                        taxRankMaps[r].put(taxon, lineage);
                    }

                    // sum intensities per taxon
                    double[] current = taxIntensityMaps[r].get(taxon);
                    if (current == null) {
                        current = new double[genomeIntensities.length];
                    }
                    for (int j = 0; j < current.length; j++) {
                        current[j] += genomeIntensities[j];
                    }
                    taxIntensityMaps[r].put(taxon, current);

                    // aggregate scores per taxon
                    Double curScore = taxScoreMaps[r].get(taxon);
                    if (isDia) {
                        if (legacyDiaProduct) {
                            taxScoreMaps[r].put(taxon, curScore == null ? genomeScore : curScore * genomeScore);
                        } else {
                            taxScoreMaps[r].put(taxon, curScore == null ? genomeScore : curScore + genomeScore);
                        }
                    } else {
                        taxScoreMaps[r].put(taxon, curScore == null ? genomeScore : curScore + genomeScore);
                    }
                }
            }

            // Header (legacy Taxa.tsv)
            StringBuilder sb = new StringBuilder();
            sb.append("Name\t");
            sb.append("Minus_log(p-value)\t");
            sb.append("Rank");
            for (String rn : rankNames) {
                sb.append("\t").append(rn);
            }
            for (String in : intensityNames) {
                sb.append("\t").append(in);
            }
            taxaWriter.println(sb);

            // Deterministic output order (sort by rank then name)
            for (int r = 0; r < gtdbRankIdx.size(); r++) {
                ArrayList<String> taxa = new ArrayList<>(taxRankMaps[r].keySet());
                Collections.sort(taxa);
                for (String taxon : taxa) {
                    String[] lineage = taxRankMaps[r].get(taxon);
                    double[] inten = taxIntensityMaps[r].get(taxon);
                    double agg = taxScoreMaps[r].getOrDefault(taxon, 0.0d);

                    StringBuilder row = new StringBuilder();
                    row.append(taxon).append("\t");
                    if (isDia && legacyDiaProduct) {
                        // legacy DIA: product(S) -> -log10(product)
                        if (agg == 0.0d) {
                            row.append("1000");
                        } else {
                            row.append(-Math.log10(agg));
                        }
                    } else {
                        // DDA legacy uses sum(-log10(p)); DIA new-style uses sum(S) if legacyDiaProduct=false
                        row.append(agg);
                    }
                    row.append("\t");
                    row.append(rankNames.get(r));

                    for (int j = 0; j < rankNames.size(); j++) {
                        if (j < lineage.length) {
                            row.append("\t").append(lineage[j] == null ? "" : lineage[j]);
                        } else {
                            row.append("\t");
                        }
                    }
                    for (int j = 0; j < inten.length; j++) {
                        row.append("\t").append(inten[j]);
                    }
                    taxaWriter.println(row);
                }
            }
        }
    }

    private static String get(String[] parts, int idx) {
        if (parts == null || idx < 0 || idx >= parts.length) return "";
        return parts[idx] == null ? "" : parts[idx];
    }

    private static int indexOf(String[] header, String colName) {
        if (header == null || colName == null) return -1;
        for (int i = 0; i < header.length; i++) {
            if (colName.equals(header[i] == null ? null : header[i].trim())) return i;
        }
        return -1;
    }

    private static double parseDoubleSafe(String[] parts, int idx, double defaultValue) {
        if (parts == null || idx < 0 || idx >= parts.length) return defaultValue;
        String s = parts[idx] == null ? "" : parts[idx].trim();
        if (s.isEmpty()) return defaultValue;
        try {
            return Double.parseDouble(s);
        } catch (NumberFormatException e) {
            return defaultValue;
        }
    }

    private static void writeDiaPeptide(String sourceDir, String targetDir) throws IOException {
        File sourceFile = new File(sourceDir);
        if (!sourceFile.exists()) {
            return;
        }
        
        PrintWriter writer = new PrintWriter(new File(targetDir,  "DIA_peptide.tsv"));
        writer.println("Method\tPeptide\tProteins");
        
        File[] files = sourceFile.listFiles();
        for (File file : files) {
            File peptideFile = new File(file, "final_peptides.tsv");
            HashMap<String, String> proteinSet = new HashMap<>();
            if(peptideFile.exists()){
                String method = file.getName();
                if(peptideFile.exists()) {
                    try (BufferedReader reader = new BufferedReader(new FileReader(peptideFile))) {
                        String line = reader.readLine();
                        String[] title = line.split("\t");
                        int peptideIndex = -1;
                        int proteinsIndex = -1;
                        for(int i = 0; i < title.length; i++) {
                            if(title[i].equals("Base sequence")) {
                                peptideIndex = i;
                            }
                            if(title[i].equals("Proteins")) {
                                proteinsIndex = i;
                            }
                        }
                        while ((line = reader.readLine()) != null) {
                            String[] parts = line.split("\t");
                            String peptide = parts[peptideIndex];
                            String proteins = parts[proteinsIndex];
                            proteinSet.put(peptide, proteins);
                        }
                    }
                }
                for(String peptide : proteinSet.keySet()) {
                    writer.println(method + "\t" + peptide + "\t" + proteinSet.get(peptide));
                }
            } else {
                System.out.println(file.getName());
            }
        }

        writer.close();
    }
    
    private static void writeDdaPeptide(String sourceDir, String targetDir) throws IOException {
        File sourceFile = new File(sourceDir);
        if (!sourceFile.exists()) {
            return;
        }
        
        PrintWriter writer = new PrintWriter(new File(targetDir,  "peptide.tsv"));
        writer.println("Method\tPeptide\tProteins");
        
        File[] files = sourceFile.listFiles();
        for (File file : files) {
            HashMap<String, String> proteinSet = new HashMap<>();
            if(file.getName().startsWith("pfind_MAG_search_")){
                String method = file.getName().substring("pfind_MAG_search_".length());
                File peptideFile = new File(file, "final_peptides.tsv");
                if(peptideFile.exists()) {
                    try (BufferedReader reader = new BufferedReader(new FileReader(peptideFile))) {
                        String line = reader.readLine();
                        String[] title = line.split("\t");
                        int peptideIndex = -1;
                        int proteinsIndex = -1;
                        for(int i = 0; i < title.length; i++) {
                            if(title[i].equals("Base sequence")) {
                                peptideIndex = i;
                            }
                            if(title[i].equals("Proteins")) {
                                proteinsIndex = i;
                            }
                        }
                        while ((line = reader.readLine()) != null) {
                            String[] parts = line.split("\t");
                            String peptide = parts[peptideIndex];
                            String proteins = parts[proteinsIndex];
                            proteinSet.put(peptide, proteins);
                        }
                    }
                }
                for(String peptide : proteinSet.keySet()) {
                    writer.println(method + "\t" + peptide + "\t" + proteinSet.get(peptide));
                }
            } else {
                System.out.println(file.getName());
            }
        }

        writer.close();
    }

    private static void writeDdaSelectedGenomes(String sourceDir, String targetDir) throws IOException {
        File sourceFile = new File(sourceDir);
        if (!sourceFile.exists()) {
            return;
        }

        String method = sourceFile.getName().substring("pfind_MAG_search_".length());
        PrintWriter writer = new PrintWriter(new File(targetDir, method + "_genome_selected.tsv"));
        writer.println("Method\tRaw File\tGenome");
        
        File[] files = sourceFile.listFiles();
        for (File file : files) {
            
            File hapFile = new File(file, "hap");
            if (hapFile.isDirectory()) {
                File selectedFile = new File(hapFile, "marker_genome_selected.tsv");
                if(!selectedFile.exists()) {
                    selectedFile = new File(hapFile, "genome_selected.tsv");
                }
                if (selectedFile.exists()) {
                    System.out.println(selectedFile.getName());
                    try (BufferedReader reader = new BufferedReader(new FileReader(selectedFile))) {
                        String line = reader.readLine();
                        while ((line = reader.readLine()) != null) {
                            if (line.startsWith("MGYG")) {
                                String[] parts = line.split("\t");
                            String genome = parts[0];
                            writer.println(method + "\t" + file.getName() + "\t" + genome);
                            }
                        }
                    }
                }
            } else {
                System.out.println(file.getName());
            }
        }

        writer.close();
    }
    
	private static void writeDiaSelectedGenomes(String sourceDir, String targetDir) throws IOException {
		File sourceFile = new File(sourceDir);
		if (!sourceFile.exists()) {
			return;
		}

		File[] files = sourceFile.listFiles();
		for (File file : files) {
			String method = file.getName();
			File seedFile = new File(file, "seed");
			if (seedFile.isDirectory() && seedFile.exists()) {
				PrintWriter writer = new PrintWriter(new File(targetDir, method + "_genome_selected.tsv"));
				writer.println("Method\tRaw File\tGenome");

				int genomeCount = 0;
				for (File roundFile : seedFile.listFiles()) {
					if (roundFile.getName().startsWith("round")) {
						genomeCount += 5;
					}
				}

				if (genomeCount == 0) {
					File deepFile = new File(seedFile, "genome_selected.tsv");
					System.out.println(deepFile.getName());
					try (BufferedReader reader = new BufferedReader(new FileReader(deepFile))) {
						String line = reader.readLine();
						while ((line = reader.readLine()) != null) {
							if (line.startsWith("MGYG")) {
								String[] parts = line.split("\t");
								String genome = parts[0];
								writer.println(method + "\t" + file.getName() + "\t" + genome);
							}
						}
					}
				} else {
					File deepFile = new File(seedFile, "genome_ranked.tsv");
					System.out.println(deepFile.getName());
					try (BufferedReader reader = new BufferedReader(new FileReader(deepFile))) {
						String line = reader.readLine();
						int count = 0;
						while ((line = reader.readLine()) != null) {
							String[] parts = line.split("\t");
							if (parts.length > 1 && parts[1].startsWith("MGYG")) {
								String genome = parts[1];
								writer.println(method + "\t" + file.getName() + "\t" + genome);
								count++;
								if (count == genomeCount) {
									break;
								}
							}
						}
					}
				}
				writer.close();
			}

			File markerFile = new File(file, "marker");
			if (markerFile.isDirectory() && markerFile.exists()) {
				PrintWriter writer = new PrintWriter(new File(targetDir, method + "_genome_selected.tsv"));
				writer.println("Method\tRaw File\tGenome");

				int genomeCount = 0;
				for (File roundFile : markerFile.listFiles()) {
					if (roundFile.getName().startsWith("round")) {
						genomeCount += 5;
					}
				}

				if (genomeCount == 0) {
					File deepFile = new File(markerFile, "genome_selected.tsv");
					System.out.println(deepFile.getName());
					try (BufferedReader reader = new BufferedReader(new FileReader(deepFile))) {
						String line = reader.readLine();
						while ((line = reader.readLine()) != null) {
							if (line.startsWith("MGYG")) {
								String[] parts = line.split("\t");
								String genome = parts[0];
								writer.println(method + "\t" + file.getName() + "\t" + genome);
							}
						}
					}
				} else {
					File deepFile = new File(markerFile, "genome_ranked.tsv");
					System.out.println(deepFile.getName());
					try (BufferedReader reader = new BufferedReader(new FileReader(deepFile))) {
						String line = reader.readLine();
						int count = 0;
						while ((line = reader.readLine()) != null) {
							String[] parts = line.split("\t");
							if (parts.length > 1 && parts[1].startsWith("MGYG")) {
								String genome = parts[1];
								writer.println(method + "\t" + file.getName() + "\t" + genome);
								count++;
								if (count == genomeCount) {
									break;
								}
							}
						}
					}
				}
				writer.close();
			}
		}
	}

    private static void writeFinalGenomes(String sourceDir, String targetDir) throws IOException {
        File sourceFile = new File(sourceDir);
        if (!sourceFile.exists()) {
            return;
        }
        
        PrintWriter writer = new PrintWriter(new File(targetDir,  "DDA_genome.tsv"));
        writer.println("Method\tGenome\tSelection_label");
        File[] files = sourceFile.listFiles();
        for (File file : files) {
            if (file.getName().startsWith("pfind_MAG_search_")) {
                String method = file.getName().substring("pfind_MAG_search_".length());
                File taxaFolder = new File(file, "taxonomy_analysis");
                if(taxaFolder.exists()) {
                    File genomeFile = new File(taxaFolder, "Genome.tsv");
                    if(genomeFile.exists()) {
                        try (BufferedReader reader = new BufferedReader(new FileReader(genomeFile))) {
                            String headerLine = reader.readLine();
                            int selIdx = 4;
                            if (headerLine != null) {
                                String[] h = headerLine.split("\t");
                                for (int c = 0; c < h.length; c++) {
                                    if ("Selection_label".equals(h[c].trim())) {
                                        selIdx = c;
                                        break;
                                    }
                                }
                            }
                            String line;
                            while ((line = reader.readLine()) != null) {
                                String[] parts = line.split("\t");
                                String genome = parts[0];
                                String sel = selIdx < parts.length ? parts[selIdx] : "";
                                writer.println(method + "\t" + genome + "\t" + sel);
                            }
                        }
                    }
                }
            } else if (file.getName().startsWith("marker") || file.getName().startsWith("seed")) {
                String method = file.getName();
                File taxaFolder = new File(file, "taxonomy_analysis");
                if(taxaFolder.exists()) {
                    File genomeFile = new File(taxaFolder, "Genome.tsv");
                    if(genomeFile.exists()) {
                        try (BufferedReader reader = new BufferedReader(new FileReader(genomeFile))) {
                            String headerLine = reader.readLine();
                            int selIdx = 4;
                            if (headerLine != null) {
                                String[] h = headerLine.split("\t");
                                for (int c = 0; c < h.length; c++) {
                                    if ("Selection_label".equals(h[c].trim())) {
                                        selIdx = c;
                                        break;
                                    }
                                }
                            }
                            String line;
                            while ((line = reader.readLine()) != null) {
                                String[] parts = line.split("\t");
                                String genome = parts[0];
                                String sel = selIdx < parts.length ? parts[selIdx] : "";
                                writer.println(method + "\t" + genome + "\t" + sel);
                            }
                        }
                    }
                }
            } else {
                System.out.println(file.getName());
            }
        }
        writer.close();
    }

    private static void combineFasta(String fasta1, String fasta2, String output) throws IOException {
        try (PrintWriter writer = new PrintWriter(output)) {
            BufferedReader reader1 = new BufferedReader(new FileReader(fasta1));
            BufferedReader reader2 = new BufferedReader(new FileReader(fasta2));
            HashMap<String, String> fastaMap = new HashMap<>();
            String line1 = null;
            while ((line1 = reader1.readLine()) != null) {
                if (line1.startsWith(">")) {
                    fastaMap.put(reader1.readLine(), line1.substring(1));
                }
            }
            String line2 = null;
            while ((line2 = reader2.readLine()) != null) {
                if (line2.startsWith(">")) {
                    fastaMap.put(reader2.readLine(), line2.substring(1));
                }
            }
            reader1.close();
            reader2.close();
            for (String key : fastaMap.keySet()) {
                writer.println(">" + fastaMap.get(key));
                writer.println(key);
            }
            writer.close();
        }
    }

    private static void combinePerFileParquet(String folder1, String folder2, String destFolder) throws IOException {
        Pattern parquetPattern = Pattern.compile("^per_file_(\\d+)\\.parquet$");

        File folder1File = new File(folder1);
        File folder2File = new File(folder2);
        File destFolderFile = new File(destFolder);

        if (!folder1File.isDirectory()) {
            throw new IOException("folder1 is not a valid directory: " + folder1);
        }
        if (!folder2File.isDirectory()) {
            throw new IOException("folder2 is not a valid directory: " + folder2);
        }
        if (!destFolderFile.exists() && !destFolderFile.mkdirs()) {
            throw new IOException("Cannot create destination directory: " + destFolder);
        }
        if (!destFolderFile.isDirectory()) {
            throw new IOException("destFolder is not a directory: " + destFolder);
        }

        List<File> folder1ParquetFiles = new ArrayList<>();
        List<File> folder2ParquetFiles = new ArrayList<>();

        File[] folder1Files = folder1File.listFiles();
        if (folder1Files != null) {
            for (File file : folder1Files) {
                if (file.isFile() && parquetPattern.matcher(file.getName()).matches()) {
                    folder1ParquetFiles.add(file);
                }
            }
        }

        File[] folder2Files = folder2File.listFiles();
        if (folder2Files != null) {
            for (File file : folder2Files) {
                if (file.isFile() && parquetPattern.matcher(file.getName()).matches()) {
                    folder2ParquetFiles.add(file);
                }
            }
        }

        folder1ParquetFiles.sort(Comparator.comparingInt(file -> extractParquetIndex(file.getName(), parquetPattern)));
        folder2ParquetFiles.sort(Comparator.comparingInt(file -> extractParquetIndex(file.getName(), parquetPattern)));

        int maxFolder1Index = -1;
        for (File file : folder1ParquetFiles) {
            int idx = extractParquetIndex(file.getName(), parquetPattern);
            Path targetPath = new File(destFolderFile, file.getName()).toPath();
            Files.copy(file.toPath(), targetPath, StandardCopyOption.REPLACE_EXISTING);
            if (idx > maxFolder1Index) {
                maxFolder1Index = idx;
            }
        }

        int nextIndex = maxFolder1Index + 1;
        for (File file : folder2ParquetFiles) {
            String newName = "per_file_" + nextIndex + ".parquet";
            Path targetPath = new File(destFolderFile, newName).toPath();
            Files.copy(file.toPath(), targetPath, StandardCopyOption.REPLACE_EXISTING);
            nextIndex++;
        }

        int expectedCount = folder1ParquetFiles.size() + folder2ParquetFiles.size();
        int actualCount = 0;
        File[] destFiles = destFolderFile.listFiles();
        if (destFiles != null) {
            for (File file : destFiles) {
                if (file.isFile() && parquetPattern.matcher(file.getName()).matches()) {
                    actualCount++;
                }
            }
        }

        if (actualCount != expectedCount) {
            throw new IOException("Parquet file count mismatch. Expected: " + expectedCount + ", actual: " + actualCount);
        }
    }

    private static int extractParquetIndex(String fileName, Pattern parquetPattern) {
        Matcher matcher = parquetPattern.matcher(fileName);
        if (!matcher.matches()) {
            throw new IllegalArgumentException("Invalid per_file parquet name: " + fileName);
        }
        return Integer.parseInt(matcher.group(1));
    }

    private static void writePeptideDetectabilityDistribution(
            String peptideListFile,
            String detectabilityFolder,
            String outputFile) throws IOException {
        File peptideFile = new File(peptideListFile);
        File detectFolder = new File(detectabilityFolder);
        File outFile = new File(outputFile);

        if (!peptideFile.isFile()) {
            throw new IOException("Peptide list file does not exist: " + peptideListFile);
        }
        if (!detectFolder.isDirectory()) {
            throw new IOException("Detectability folder is not valid: " + detectabilityFolder);
        }
        File outParent = outFile.getParentFile();
        if (outParent != null && !outParent.exists() && !outParent.mkdirs()) {
            throw new IOException("Cannot create output folder: " + outParent.getAbsolutePath());
        }

        List<String> peptideOrder = new ArrayList<>();
        Set<String> peptideSet = new HashSet<>();
        HashMap<String, Double> peptideMaxDetectability = new HashMap<>();

        try (BufferedReader reader = new BufferedReader(new FileReader(peptideFile))) {
            String header = reader.readLine();
            if (header == null) {
                throw new IOException("Peptide list file is empty: " + peptideListFile);
            }

            String[] headerParts = header.split("\t", -1);
            int peptideIndex = -1;
            for (int i = 0; i < headerParts.length; i++) {
                if ("Stripped.Sequence".equals(headerParts[i].trim())) {
                    peptideIndex = i;
                    break;
                }
            }
            if (peptideIndex < 0) {
                throw new IOException("Column 'Stripped.Sequence' not found in: " + peptideListFile);
            }

            String line;
            while ((line = reader.readLine()) != null) {
                if (line.trim().isEmpty()) {
                    continue;
                }
                String[] parts = line.split("\t", -1);
                if (peptideIndex >= parts.length) {
                    continue;
                }
                String peptide = parts[peptideIndex].trim();
                if (peptide.isEmpty()) {
                    continue;
                }
                if (peptideSet.add(peptide)) {
                    peptideOrder.add(peptide);
                }
            }
        }

        try (Stream<Path> fileStream = Files.walk(detectFolder.toPath())) {
            fileStream.filter(Files::isRegularFile).forEach(filePath -> {
                try (BufferedReader reader = new BufferedReader(new FileReader(filePath.toFile()))) {
                    String header = reader.readLine();
                    if (header == null) {
                        return;
                    }

                    String[] headers = header.split("\t", -1);
                    int peptideIndex = -1;
                    int detectabilityIndex = -1;
                    for (int i = 0; i < headers.length; i++) {
                        String column = headers[i].trim();
                        if ("Peptide sequence".equals(column)) {
                            peptideIndex = i;
                        } else if ("Peptide detectability".equals(column)) {
                            detectabilityIndex = i;
                        }
                    }
                    if (peptideIndex < 0 || detectabilityIndex < 0) {
                        return;
                    }

                    String line;
                    while ((line = reader.readLine()) != null) {
                        if (line.trim().isEmpty()) {
                            continue;
                        }
                        String[] parts = line.split("\t", -1);
                        if (peptideIndex >= parts.length || detectabilityIndex >= parts.length) {
                            continue;
                        }

                        String peptide = parts[peptideIndex].trim();
                        if (!peptideSet.contains(peptide)) {
                            continue;
                        }

                        String detectabilityText = parts[detectabilityIndex].trim();
                        if (detectabilityText.isEmpty()) {
                            continue;
                        }

                        try {
                            double detectability = Double.parseDouble(detectabilityText);
                            Double currentMax = peptideMaxDetectability.get(peptide);
                            if (currentMax == null || detectability > currentMax) {
                                peptideMaxDetectability.put(peptide, detectability);
                            }
                        } catch (NumberFormatException ignored) {
                            // Skip malformed detectability values.
                        }
                    }
                } catch (IOException ignored) {
                    // Skip unreadable files and continue scanning.
                }
            });
        }

        System.out.println("Total peptides in list: " + peptideOrder.size());
        System.out.println("Peptides matched with detectability: " + peptideMaxDetectability.size());
        for (int thresholdStep = 1; thresholdStep <= 19; thresholdStep++) {
            double threshold = thresholdStep * 0.05d;
            int count = 0;
            for (String peptide : peptideOrder) {
                double detectability = peptideMaxDetectability.getOrDefault(peptide, 0.0d);
                if (detectability > threshold) {
                    count++;
                }
            }
            System.out.println(String.format("Accumulated peptide number (>%1$.2f): %2$d", threshold, count));
        }

        try (PrintWriter writer = new PrintWriter(outFile)) {
            writer.println("Peptide\tPeptide detectability");
            for (String peptide : peptideOrder) {
                double detectability = peptideMaxDetectability.getOrDefault(peptide, 0.0d);
                writer.println(peptide + "\t" + detectability);
            }
        }
    }

    private static void getNoneZeroIntensityPeptides(String intensityFile) throws IOException {
        try (BufferedReader reader = new BufferedReader(new FileReader(intensityFile))) {
            String header = reader.readLine();
            if (header == null) {
                throw new IOException("Intensity file is empty: " + intensityFile);
            }

            String[] headerParts = header.split("\t", -1);
            ArrayList<Integer> intensityIndices = new ArrayList<>();
            ArrayList<String> intensityNames = new ArrayList<>();
            for (int i = 0; i < headerParts.length; i++) {
//                if (headerParts[i].trim().startsWith("Intensity ")) {
//                    intensityIndices.add(i);
//                    intensityNames.add(headerParts[i].substring("Intensity ".length()));
//                }
                if (headerParts[i].trim().startsWith("Z:")) {
                    intensityIndices.add(i);
                    String[] intensityNameParts = headerParts[i].split("\\\\", -1);
                    String intensityName = intensityNameParts[intensityNameParts.length - 1];
                    intensityNames.add(intensityName);
                }
            }

            int[] intensityCounts = new int[intensityNames.size()];
            String line = null;
            while ((line = reader.readLine()) != null) {
                if (line.trim().isEmpty()) {
                    continue;
                }
                String[] parts = line.split("\t", -1);
                for (int i = 0; i < intensityIndices.size(); i++) {
                    if(parts[intensityIndices.get(i)].trim().isEmpty()) {
                        continue;
                    }
                    double intensity = Double.parseDouble(parts[intensityIndices.get(i)]);
                    if (intensity > 0) {
                        intensityCounts[i]++;
                    }
                }
            }
            for (int i = 0; i < intensityNames.size(); i++) {
                System.out.println(intensityNames.get(i) + "\t" + intensityCounts[i]);
            }
            reader.close();
        }
    }

    private static void comparePeptideLevel3(String file1, String file2, String file3) throws IOException {
        HashSet<String> peptideSet1 = new HashSet<>();
        HashSet<String> peptideSet2 = new HashSet<>();
        HashSet<String> peptideSet3 = new HashSet<>();

        // File1 only: genome distribution for peptides from my result (protein column idx=5).
        // genome -> unique peptide set (peptides can map to multiple genomes if proteins list contains multiple genomes)
        HashMap<String, HashSet<String>> genomeToPeptides1 = new HashMap<>();
        // Group index: 0,1,2,3 (3 means >=3 missed cleavages)
        HashMap<Integer, HashSet<String>> peptideGroups1 = new HashMap<>();
        HashMap<Integer, HashSet<String>> peptideGroups2 = new HashMap<>();
        HashMap<Integer, HashSet<String>> peptideGroups3 = new HashMap<>();
        for (int i = 0; i <= 3; i++) {
            peptideGroups1.put(i, new HashSet<>());
            peptideGroups2.put(i, new HashSet<>());
            peptideGroups3.put(i, new HashSet<>());
        }

        try (BufferedReader reader1 = new BufferedReader(new FileReader(file1))) {
            String header1 = reader1.readLine();
            if (header1 == null) {
                throw new IOException("File 1 is empty: " + file1);
            }

            String line1 = null;
            while ((line1 = reader1.readLine()) != null) {
                if (line1.trim().isEmpty()) {
                    continue;
                }
                String[] parts1 = line1.split("\t", -1);
                if (parts1.length <= 1) {
                    continue;
                }
                String peptide1 = parts1[1] == null ? "" : parts1[1].trim();
                if (peptide1.isEmpty()) {
                    continue;
                }
                String pep1 = normalizeIL(peptide1);
                peptideSet1.add(pep1);
                int group1 = toMissedCleavageGroup(countMissedCleavagesKR(pep1));
                peptideGroups1.get(group1).add(pep1);

                // Protein column is the 6th column (idx=5): use it to infer genome(s).
                // Expected protein IDs like "MGYGxxxx_protein" (possibly multiple proteins separated by ';' or ',').
                if (parts1.length > 5) {
                    String proteinsCell = parts1[5] == null ? "" : parts1[5].trim();
                    if (!proteinsCell.isEmpty()) {
                        // Split common multi-protein separators.
                        String[] prots = proteinsCell.split("[;,]");
                        for (String p : prots) {
                            if (p == null) continue;
                            String pro = p.trim();
                            if (pro.isEmpty()) continue;

                            // Some exports include multiple tokens separated by whitespace; keep only first token.
                            int sp = pro.indexOf(' ');
                            if (sp > 0) pro = pro.substring(0, sp);

                            if (pro.startsWith("REV_")) pro = pro.substring(4);

                            // Genome id is prefix before first '_' (matches MarkerGenomeSelector.extractGenomeId behavior).
                            int us = pro.indexOf('_');
                            if (us <= 0) continue;
                            String genome = pro.substring(0, us).trim();
                            if (genome.isEmpty()) continue;

                            genomeToPeptides1.computeIfAbsent(genome, k -> new HashSet<>()).add(pep1);
                        }
                    }
                }
            }
        }
    
        try (BufferedReader reader2 = new BufferedReader(new FileReader(file2))) {
            String header2 = reader2.readLine();
            if (header2 == null) {
                throw new IOException("File 2 is empty: " + file2);
            }
            String[] h2 = header2.split("\t", -1);
            int peptideCol2 = detectPeptideColumnIndex(h2, 6);

            String line2 = null;
            while ((line2 = reader2.readLine()) != null) {
                if (line2.trim().isEmpty()) {
                    continue;
                }
                String[] parts2 = line2.split("\t", -1);
                if (peptideCol2 < 0 || parts2.length <= peptideCol2) {
                    continue;
                }
                String peptide2 = parts2[peptideCol2] == null ? "" : parts2[peptideCol2].trim();
                if (peptide2.isEmpty()) {
                    continue;
                }
                String pep2 = normalizeIL(peptide2);
                peptideSet2.add(pep2);
                int group2 = toMissedCleavageGroup(countMissedCleavagesKR(pep2));
                peptideGroups2.get(group2).add(pep2);
            }
        }

        try (BufferedReader reader3 = new BufferedReader(new FileReader(file3))) {
            String header3 = reader3.readLine();
            if (header3 == null) {
                throw new IOException("File 3 is empty: " + file3);
            }
            String[] h3 = header3.split("\t", -1);
            int peptideCol3 = detectPeptideColumnIndex(h3, 5);

            String line3 = null;
            while ((line3 = reader3.readLine()) != null) {
                if (line3.trim().isEmpty()) {
                    continue;
                }
                String[] parts3 = line3.split("\t", -1);
                if (peptideCol3 < 0 || parts3.length <= peptideCol3) {
                    continue;
                }
                String peptide3 = parts3[peptideCol3] == null ? "" : parts3[peptideCol3].trim();
                if (peptide3.isEmpty()) {
                    continue;
                }
                String pep3 = normalizeIL(peptide3);
                peptideSet3.add(pep3);
                int group3 = toMissedCleavageGroup(countMissedCleavagesKR(pep3));
                peptideGroups3.get(group3).add(pep3);
            }
        }

        System.out.println("Peptides in file 1: " + peptideSet1.size());
        System.out.println("Peptides in file 2: " + peptideSet2.size());
        System.out.println("Peptides in file 3: " + peptideSet3.size());

        // Export peptides identified by paper but missing from my result:
        // peptides present in file2 OR file3, but NOT in file1.
        {
            HashSet<String> paperPeptides = new HashSet<>(peptideSet2);
            paperPeptides.addAll(peptideSet3);
            paperPeptides.removeAll(peptideSet1);

            ArrayList<String> missing = new ArrayList<>(paperPeptides);
            Collections.sort(missing);

            File outDir = new File(file1).getParentFile();
            if (outDir == null) outDir = new File(".");
            File outFile = new File(outDir, "peptides_in_paper_missing_from_my_result.tsv");

            try (PrintWriter w = new PrintWriter(outFile)) {
                w.println("Peptide\tInFile2\tInFile3");
                for (String pep : missing) {
                    w.println(pep + "\t" + (peptideSet2.contains(pep) ? "1" : "0") + "\t" + (peptideSet3.contains(pep) ? "1" : "0"));
                }
            }
            System.out.println("Exported missing peptides: " + missing.size() + "\t" + outFile.getAbsolutePath());
        }

        // Print genome distribution for peptides from file1 (unique peptides per genome).
        if (!genomeToPeptides1.isEmpty()) {
            ArrayList<Map.Entry<String, HashSet<String>>> entries = new ArrayList<>(genomeToPeptides1.entrySet());
            entries.sort((a, b) -> {
                int ca = a.getValue() == null ? 0 : a.getValue().size();
                int cb = b.getValue() == null ? 0 : b.getValue().size();
                if (cb != ca) return Integer.compare(cb, ca);
                return a.getKey().compareTo(b.getKey());
            });
            System.out.println("---- File1 genome distribution (unique peptides; protein col idx=5) ----");
            for (Map.Entry<String, HashSet<String>> e : entries) {
                System.out.println(e.getKey() + "\t" + (e.getValue() == null ? 0 : e.getValue().size()));
            }
        } else {
            System.out.println("---- File1 genome distribution ----");
            System.out.println("No genome/protein info parsed from file1 column 6 (idx=5).");
        }

        System.out.println("---- Venn3 counts (overall) ----");
        printVenn3Counts(peptideSet1, peptideSet2, peptideSet3);

        System.out.println("---- Venn3 counts by missed cleavage group (KR, excluding KP/RP) ----");
        for (int g = 0; g <= 3; g++) {
            String label = (g <= 2) ? String.valueOf(g) : ">=3";
            System.out.println("MissedCleavage=" + label);
            printVenn3Counts(peptideGroups1.get(g), peptideGroups2.get(g), peptideGroups3.get(g));
        }
    }

    /**
     * Compare 4 peptide files and export UpSet plot counts (intersection sizes).
     * Peptides are normalized by equating I/L (I -> L).
     *
     * File formats:
     * - file1: peptide column fixed at idx=1 (as in comparePeptideLevel3).
     * - file2: peptide column detected (fallback idx=6).
     * - file3: peptide column detected (fallback idx=5).
     * - file4: same format as file2 (detected, fallback idx=6).
     *
     * Output: a TSV next to file1 named "upset4_peptide_counts.tsv".
     */
    private static void comparePeptideLevel4Upset(String file1, String file2, String file3, String file4) throws IOException {
        HashSet<String> s1 = new HashSet<>();
        HashSet<String> s2 = new HashSet<>();
        HashSet<String> s3 = new HashSet<>();
        HashSet<String> s4 = new HashSet<>();

        // file1 (my result): peptide column idx=1
        try (BufferedReader r1 = new BufferedReader(new FileReader(file1))) {
            String h = r1.readLine();
            if (h == null) throw new IOException("File 1 is empty: " + file1);
            String line;
            while ((line = r1.readLine()) != null) {
                if (line.trim().isEmpty()) continue;
                String[] parts = line.split("\t", -1);
                if (parts.length <= 1) continue;
                String pep = parts[1] == null ? "" : parts[1].trim();
                if (pep.isEmpty()) continue;
                s1.add(normalizeIL(pep));
            }
        }
        System.out.println("s1\t"+s1.size());
        // file2
        try (BufferedReader r2 = new BufferedReader(new FileReader(file2))) {
            String header = r2.readLine();
            if (header == null) throw new IOException("File 2 is empty: " + file2);
            int pepCol = detectPeptideColumnIndex(header.split("\t", -1), 6);
            String line;
            while ((line = r2.readLine()) != null) {
                if (line.trim().isEmpty()) continue;
                String[] parts = line.split("\t", -1);
                if (pepCol < 0 || parts.length <= pepCol) continue;
                String pep = parts[pepCol] == null ? "" : parts[pepCol].trim();
                if (pep.isEmpty()) continue;
                s2.add(normalizeIL(pep));
            }
        }
        System.out.println("s2\t"+s2.size());
        // file3
        try (BufferedReader r3 = new BufferedReader(new FileReader(file3))) {
            String header = r3.readLine();
            if (header == null) throw new IOException("File 3 is empty: " + file3);
            int pepCol = detectPeptideColumnIndex(header.split("\t", -1), 5);
            String line;
            while ((line = r3.readLine()) != null) {
                if (line.trim().isEmpty()) continue;
                String[] parts = line.split("\t", -1);
                if (pepCol < 0 || parts.length <= pepCol) continue;
                String pep = parts[pepCol] == null ? "" : parts[pepCol].trim();
                if (pep.isEmpty()) continue;
                s3.add(normalizeIL(pep));
            }
        }
        System.out.println("s3\t"+s3.size());
        // file4 (same as file2)
        try (BufferedReader r4 = new BufferedReader(new FileReader(file4))) {
            String header = r4.readLine();
            if (header == null) throw new IOException("File 4 is empty: " + file4);
            int pepCol = detectPeptideColumnIndex(header.split("\t", -1), 6);
            String line;
            while ((line = r4.readLine()) != null) {
                if (line.trim().isEmpty()) continue;
                String[] parts = line.split("\t", -1);
                if (pepCol < 0 || parts.length <= pepCol) continue;
                String pep = parts[pepCol] == null ? "" : parts[pepCol].trim();
                if (pep.isEmpty()) continue;
                s4.add(normalizeIL(pep));
            }
        }
        System.out.println("s4\t"+s4.size());
        // Union
        HashSet<String> union = new HashSet<>(s1);
        union.addAll(s2);
        union.addAll(s3);
        union.addAll(s4);
        System.out.println("union\t"+union.size());
        // Export peptide list: only in file1 AND file4, and NOT in file2 and file3.
        // (I/L equated already via normalizeIL during read)
        HashSet<String> only14 = new HashSet<>(s1);
        only14.retainAll(s4);
        only14.removeAll(s2);
        only14.removeAll(s3);
        System.out.println("only14\t"+only14.size());
        // mask: bit0=file1, bit1=file2, bit2=file3, bit3=file4
        HashMap<Integer, Integer> maskCounts = new HashMap<>();
        for (String pep : union) {
            int m = 0;
            if (s1.contains(pep)) m |= 1;
            if (s2.contains(pep)) m |= 2;
            if (s3.contains(pep)) m |= 4;
            if (s4.contains(pep)) m |= 8;
            if (m == 0) continue;
            maskCounts.put(m, maskCounts.getOrDefault(m, 0) + 1);
        }

        ArrayList<Integer> masks = new ArrayList<>(maskCounts.keySet());
        masks.sort((a, b) -> {
            int ca = maskCounts.getOrDefault(a, 0);
            int cb = maskCounts.getOrDefault(b, 0);
            if (cb != ca) return Integer.compare(cb, ca);
            return Integer.compare(a, b);
        });

        File outDir = new File(file1).getParentFile();
        if (outDir == null) outDir = new File(".");
        File outFile = new File(outDir, "upset4_peptide_counts.tsv");
        File outOnly14 = new File(outDir, "peptides_only_file1_and_file4.tsv");
        File outUnion = new File(outDir, "peptides_union_4files.tsv");
        File outCommon = new File(outDir, "peptides_common_4files.tsv");

        try (PrintWriter w = new PrintWriter(outFile)) {
            w.println("Mask\tInFile1\tInFile2\tInFile3\tInFile4\tCount");
            for (int m : masks) {
                int c = maskCounts.getOrDefault(m, 0);
                w.println(m + "\t" + ((m & 1) != 0 ? "1" : "0")
                        + "\t" + ((m & 2) != 0 ? "1" : "0")
                        + "\t" + ((m & 4) != 0 ? "1" : "0")
                        + "\t" + ((m & 8) != 0 ? "1" : "0")
                        + "\t" + c);
            }
        }

        ArrayList<String> only14List = new ArrayList<>(only14);
        Collections.sort(only14List);
        try (PrintWriter w = new PrintWriter(outOnly14)) {
            w.println("Peptide");
            for (String pep : only14List) {
                w.println(pep);
            }
        }

        ArrayList<String> unionList = new ArrayList<>(union);
        Collections.sort(unionList);
        try (PrintWriter w = new PrintWriter(outUnion)) {
            w.println("Peptide\tInFile1\tInFile2\tInFile3\tInFile4");
            for (String pep : unionList) {
                w.println(pep
                        + "\t" + (s1.contains(pep) ? "1" : "0")
                        + "\t" + (s2.contains(pep) ? "1" : "0")
                        + "\t" + (s3.contains(pep) ? "1" : "0")
                        + "\t" + (s4.contains(pep) ? "1" : "0"));
            }
        }

        HashSet<String> common = new HashSet<>(s1);
        common.retainAll(s2);
        common.retainAll(s3);
        common.retainAll(s4);
        ArrayList<String> commonList = new ArrayList<>(common);
        Collections.sort(commonList);
        try (PrintWriter w = new PrintWriter(outCommon)) {
            w.println("Peptide");
            for (String pep : commonList) {
                w.println(pep);
            }
        }
        System.out.println("common\t"+common.size());
        System.out.println("---- UpSet4 summary (I/L equated) ----");
        System.out.println("File1 peptides: " + s1.size());
        System.out.println("File2 peptides: " + s2.size());
        System.out.println("File3 peptides: " + s3.size());
        System.out.println("File4 peptides: " + s4.size());
        System.out.println("Union peptides: " + union.size());
        System.out.println("Exported UpSet counts: " + masks.size() + "\t" + outFile.getAbsolutePath());
        System.out.println("Exported peptides only in file1&file4 (not in file2/file3): " + only14List.size()
                + "\t" + outOnly14.getAbsolutePath());
        System.out.println("Exported union peptide list (4 files): " + unionList.size()
                + "\t" + outUnion.getAbsolutePath());
        System.out.println("Exported common peptide list (all 4 files): " + commonList.size()
                + "\t" + outCommon.getAbsolutePath());
    }

    private static void printVenn3Counts(Set<String> a, Set<String> b, Set<String> c) {
        HashSet<String> aOnly = new HashSet<>(a);
        aOnly.removeAll(b);
        aOnly.removeAll(c);

        HashSet<String> bOnly = new HashSet<>(b);
        bOnly.removeAll(a);
        bOnly.removeAll(c);

        HashSet<String> cOnly = new HashSet<>(c);
        cOnly.removeAll(a);
        cOnly.removeAll(b);

        HashSet<String> ab = new HashSet<>(a);
        ab.retainAll(b);
        HashSet<String> abc = new HashSet<>(ab);
        abc.retainAll(c);

        HashSet<String> abOnly = new HashSet<>(ab);
        abOnly.removeAll(c);

        HashSet<String> ac = new HashSet<>(a);
        ac.retainAll(c);
        HashSet<String> acOnly = new HashSet<>(ac);
        acOnly.removeAll(b);

        HashSet<String> bc = new HashSet<>(b);
        bc.retainAll(c);
        HashSet<String> bcOnly = new HashSet<>(bc);
        bcOnly.removeAll(a);

        System.out.println("A_only\t" + aOnly.size());
        System.out.println("B_only\t" + bOnly.size());
        System.out.println("C_only\t" + cOnly.size());
        System.out.println("A∩B_only\t" + abOnly.size());
        System.out.println("A∩C_only\t" + acOnly.size());
        System.out.println("B∩C_only\t" + bcOnly.size());
        System.out.println("A∩B∩C\t" + abc.size());
    }

    private static int detectPeptideColumnIndex(String[] header, int fallbackIdx) {
        if (header == null) return fallbackIdx;
        // Common peptide column names in different exports
        String[] candidates = new String[] {
                "Base sequence",
                "Stripped.Sequence",
                "Peptide",
                "Peptide sequence",
                "Sequence", "PeptideSequence"
        };

        for (String cand : candidates) {
            for (int i = 0; i < header.length; i++) {
                if (cand.equals(header[i] == null ? null : header[i].trim())) {
                    return i;
                }
            }
        }
        return fallbackIdx;
    }

    private static int toMissedCleavageGroup(int missedCleavages) {
        if (missedCleavages <= 0) return 0;
        if (missedCleavages == 1) return 1;
        if (missedCleavages == 2) return 2;
        return 3;
    }

    /**
     * Counts missed cleavage sites for K/R cleavage with proline exception (no cut at KP/RP).
     * For a peptide sequence, this is approximated as the number of internal cleavage sites:
     * occurrences of K or R that are not the last residue and not followed by P.
     *
     * This method tolerates simple modified strings by stripping non-letter characters first.
     */
    private static int countMissedCleavagesKR(String peptide) {
        if (peptide == null) return 0;

        StringBuilder seq = new StringBuilder(peptide.length());
        for (int i = 0; i < peptide.length(); i++) {
            char c = peptide.charAt(i);
            if (Character.isLetter(c)) {
                seq.append(Character.toUpperCase(c));
            }
        }

        int len = seq.length();
        if (len <= 1) return 0;

        int missed = 0;
        for (int i = 0; i < len - 1; i++) {
            char c = seq.charAt(i);
            if (c == 'K' || c == 'R') {
                char next = seq.charAt(i + 1);
                if (next != 'P') {
                    missed++;
                }
            }
        }
        return missed;
    }
/**
 * @d
 * @param args
 * @throws IOException
 */
    public static void main(String[] args) throws IOException {

//    	comparePeptideLevel3("Z:\\PRIDE_dataset\\PXD051104_metaExpertPro\\DIA_combined\\MetaPilot\\marker_deep\\final_peptides.tsv", 
//        "Z:\\PRIDE_dataset\\PXD051104_metaExpertPro\\O480_60min_DIA_IGC_humanswiss.pr_matrix.tsv",
//        "Z:\\PRIDE_dataset\\PXD051104_metaExpertPro"
//    			+ "\\O480_60min_IGC_humanswiss_library.tsv");
    	
    	comparePeptideLevel4Upset("Z:\\PRIDE_dataset\\PXD051104_metaExpertPro\\DIA_combined\\MetaPilot\\marker_deep\\final_peptides.tsv", 
    	        "Z:\\PRIDE_dataset\\PXD051104_metaExpertPro\\O480_60min_DIA_IGC_humanswiss.pr_matrix.tsv",
    	        "Z:\\PRIDE_dataset\\PXD051104_metaExpertPro"
    	    			+ "\\O480_60min_IGC_humanswiss_library.tsv", 
    	    			"Z:\\PRIDE_dataset\\PXD051104_metaExpertPro\\DIA_combined\\MGYG000002293\\MGYG000002293.pr_matrix.tsv");
    	
//        getNoneZeroIntensityPeptides("Z:\\PRIDE_dataset\\PXD051104_metaExpertPro\\DIA_combined\\MetaPilot\\marker_deep\\host\\host.pr_matrix.tsv");
//    	exportLegacyTaxaFromGenomeTsv(new File("Z:\\PRIDE_dataset\\PXD051104_metaExpertPro\\DIA_combined\\MetaPilot\\marker_deep\\taxonomy_analysis"
//    			+ "\\Genome.tsv"), new File("Z:\\PRIDE_dataset\\PXD051104_metaExpertPro\\DIA_combined\\MetaPilot\\marker_deep\\taxonomy_analysis"
//    					+ "\\taxa_new.tsv"), false);
    	
//    	writePeptideDetectabilityDistribution("Z:\\PRIDE_dataset\\PXD051104_metaExpertPro\\DIA_combined\\temp\\combined.pr_matrix.tsv", 
//    			"C:\\Users\\daz25hak\\MetaProAI\\database\\MAG\\human-gut\\v2.0.2\\predicted", 
//    			"Z:\\PRIDE_dataset\\PXD051104_metaExpertPro\\DIA_combined\\temp\\peptide_dect.tsv");
    	
//    	writePeptideDetectabilityDistribution("Z:\\PRIDE_dataset\\PXD051792_uMetaP\\20241012_colonic_content\\MetaPilot\\marker_deep\\combined\\combined.pr_matrix.tsv", 
//   			"C:\\Users\\daz25hak\\MetaProAI\\database\\MAG\\human-gut\\v2.0.2\\predicted", 
//    			"Z:\\PRIDE_dataset\\PXD051104_metaExpertPro\\DIA_combined\\temp\\peptide_dect.tsv");
    	
//    	combinePerFileParquet("Z:\\PRIDE_dataset\\PXD051104_metaExpertPro\\GNHSF_DIA\\MetaPilot\\marker_deep\\marker", "Z:\\PRIDE_dataset\\PXD051104_metaExpertPro\\MetaPilot_DIA_MPW\\MetaPilot\\marker_deep\\marker", 
//    			"Z:\\PRIDE_dataset\\PXD051104_metaExpertPro\\DIA_combined\\MetaPilot\\marker_deep\\marker");
//        combineFasta("Z:\\PRIDE_dataset\\PXD051104_metaExpertPro\\MetaPilot_DIA_MPW\\MetaPilot\\marker_deep\\combined\\combined.fasta", 
//        "Z:\\PRIDE_dataset\\PXD051104_metaExpertPro\\GNHSF_DIA\\MetaPilot\\marker_deep\\combined\\combined.fasta", 
//        "Z:\\PRIDE_dataset\\PXD051104_metaExpertPro\\DIA_combined\\combined.fasta");
//        writeFinalGenomes("Z:\\PRIDE_dataset\\IPX0003851000\\IPX0003851004\\DIA\\MetaPilot",
//                "Z:\\PRIDE_dataset\\PXD008738_Elo_DIA\\ai_analysis\\data\\Qiao-fecal");
//                writeFinalGenomes("Z:\\PRIDE_dataset\\PXD008738_Elo_DIA\\12mix\\DDA\\MetaPilot",
//                "Z:\\PRIDE_dataset\\PXD008738_Elo_DIA\\ai_analysis\\data\\Elo-12mix");
//        writeDiaSelectedGenomes("Z:\\PRIDE_dataset\\PXD008738_Elo_DIA\\12mix\\DIA\\MetaPilot",
//                "Z:\\PRIDE_dataset\\PXD008738_Elo_DIA\\ai_analysis\\data\\Elo-12mix");
//        writeDdaSelectedGenomes("Z:\\PRIDE_dataset\\IPX0003851000\\IPX0003851001\\DDA\\MetaPilot\\"
//        		+ "pfind_MAG_search_marker2_comp", 
//        "Z:\\PRIDE_dataset\\PXD008738_Elo_DIA\\ai_analysis\\data\\Qiao-12mix");
 
//writeDiaPeptide("Z:\\PRIDE_dataset\\IPX0003851000\\IPX0003851001\\DIA\\MetaPilot",
//				"Z:\\PRIDE_dataset\\PXD008738_Elo_DIA\\ai_analysis\\data\\Qiao-12mix");
//               writeDiaPeptide("Z:\\PRIDE_dataset\\IPX0003851000\\IPX0003851004\\DIA\\MetaPilot",
//				"Z:\\PRIDE_dataset\\PXD008738_Elo_DIA\\ai_analysis\\data\\Qiao-fecal");
//                writeDiaPeptide("Z:\\PRIDE_dataset\\PXD008738_Elo_DIA\\12mix\\DIA\\MetaPilot",
//				"Z:\\PRIDE_dataset\\PXD008738_Elo_DIA\\ai_analysis\\data\\Elo-12mix");
//                writeDiaPeptide("Z:\\PRIDE_dataset\\PXD008738_Elo_DIA\\pool\\DIA\\MetaPilot",
//				"Z:\\PRIDE_dataset\\PXD008738_Elo_DIA\\ai_analysis\\data\\Elo-fecal");
/*    	
    	File[] files = (new File("Z:\\PRIDE_dataset\\PXD007819_xu_IBD\\MetaPilot\\pfind_MAG_search")).listFiles();
    	HashMap<String, Integer> map = new HashMap<String, Integer>();
    	for(File f: files) {
    		if(f.getName().startsWith("Xu_")) {
    			map.put(f.getName(), 0);
    		}
    	}
  
    	BufferedReader reader = new BufferedReader(new FileReader("Z:\\PRIDE_dataset\\PXD007819_xu_IBD\\MetaPilot\\pfind_MAG_search\\combined"
    			+ "\\psms.tsv"));
    	String line = null;
    	int[] count = new int[4];
    	HashSet<String> set = new HashSet<String>();
    	while((line=reader.readLine())!=null) {
    		String[] cs = line.split("\t");
    		if(cs[3].contains("REV_")) {
    			count[1]++;
    		}else {
    			count[0]++;
    		}
    		if(cs[3].contains("sp|")) {
    			count[2]++;
    		}else {
    			count[3]++;
    		}
    		set.add(cs[0]);
    		if(map.containsKey(cs[0])) {
    			map.put(cs[0], map.get(cs[0])+1);
    		}
    	}
            reader.close();  
        System.out.println(Arrays.toString(count)+"\t"+set.size());
        for(String filename: map.keySet()) {
        	System.out.println(filename+"\t"+map.get(filename));
        }
*/        
    }
}
