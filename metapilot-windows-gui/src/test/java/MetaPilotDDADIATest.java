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
import java.util.Set;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.stream.Stream;

import cern.colt.Arrays;

public class MetaPilotDDADIATest {

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
/**
 * @d
 * @param args
 * @throws IOException
 */
    public static void main(String[] args) throws IOException {
    	
        getNoneZeroIntensityPeptides("Z:\\PRIDE_dataset\\PXD051104_metaExpertPro\\DIA_combined\\MetaPilot\\marker_deep\\host\\host.pr_matrix.tsv");
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
