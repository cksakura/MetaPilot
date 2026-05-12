package uk.ac.quadram.metapilot.core.runner;

import uk.ac.quadram.metapilot.core.dia.DiaScoring;
import uk.ac.quadram.metapilot.core.prodb.FastaManager;
import uk.ac.quadram.metapilot.core.tools.FormatTool;
import uk.ac.quadram.metapilot.core.workflow.ProgressCallback;
import uk.ac.quadram.metapilot.dbSearch.diann.DiaNNParquetReader;
import uk.ac.quadram.metapilot.dbSearch.diann.DiaNNPrecursor;
import uk.ac.quadram.metapilot.dbSearch.diann.DiannPrMatrixReader;
import uk.ac.quadram.metapilot.quant.flashLFQ.FlashLfqQuanPeptide;
import uk.ac.quadram.metapilot.task.io.pep.MetaPeptide;
import uk.ac.quadram.metapilot.task.io.pro.MetaProtein;
import uk.ac.quadram.metapilot.task.mag.MagDbItem;
import uk.ac.quadram.metapilot.task.mag.MarkerGenomeSelector;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.math.BigInteger;
import java.nio.file.Files;
import java.nio.file.StandardCopyOption;
import java.text.DecimalFormat;
import java.util.*;

/**
 * Headless post-analysis for DIA Marker workflow (no HTML).
 *
 * Writes under resultDir:
 * - id.tsv
 * - final_peptides.tsv
 * - final_proteins.tsv
 * - all_proteins.tsv
 * - final_summary.tsv
 * - taxonomy_analysis/Genome.tsv
 * - taxonomy_analysis/Taxa.tsv
 * - taxonomy_analysis/Taxa_report.tsv
 * - taxonomy_analysis/genome_ranked.tsv (copy of combined ranking used for Selection_label)
 * - functional_annotation/* (when func_def.db is configured)
 *
 * DIA scoring is delegated to {@link DiaScoring} (shared with the Windows GUI).
 */
public final class DiaMarkerPostAnalysisTsvExporter {

    private DiaMarkerPostAnalysisTsvExporter() {
    }

    private static final int IDX_S = DiaScoring.IDX_S;
    private static final int IDX_P_EMP = DiaScoring.IDX_P_EMP;
    private static final int IDX_Q_VALUE = DiaScoring.IDX_Q_VALUE;

    public static final class Result {
        public final File idTsv;
        public final File finalPeptidesTsv;
        public final File finalProteinsTsv;
        public final File taxonomyDir;
        public final File genomeTsv;
        public final File taxaReportTsv;
        public final HashMap<String, double[]> genomeIntensityMap;
        public final HashMap<String, double[]> genomeDiaStats;

        Result(File idTsv,
               File finalPeptidesTsv,
               File finalProteinsTsv,
               File taxonomyDir,
               File genomeTsv,
               File taxaReportTsv,
               HashMap<String, double[]> genomeIntensityMap,
               HashMap<String, double[]> genomeDiaStats) {
            this.idTsv = idTsv;
            this.finalPeptidesTsv = finalPeptidesTsv;
            this.finalProteinsTsv = finalProteinsTsv;
            this.taxonomyDir = taxonomyDir;
            this.genomeTsv = genomeTsv;
            this.taxaReportTsv = taxaReportTsv;
            this.genomeIntensityMap = genomeIntensityMap;
            this.genomeDiaStats = genomeDiaStats;
        }
    }

    public static Result export(
            File resultDir,
            MagDbItem magDb,
            File combinedParquet,
            File combinedPrMatrix,
            File markerGenomeRankedTsv,
            File optionalGtdbNcbiMappingTsv,
            String markerMode,
            File funcDbFile,
            String microDbCatalogueId,
            String quantitativeMode,
            ProgressCallback cb
    ) throws Exception {
        if (cb == null) cb = new ProgressCallback.NoOp();
        if (resultDir == null) throw new IOException("resultDir is null");
        if (magDb == null) throw new IOException("magDb is null");
        if (combinedParquet == null || !combinedParquet.exists())
            throw new IOException("combined.parquet missing: " + (combinedParquet == null ? "null" : combinedParquet.getAbsolutePath()));
        if (combinedPrMatrix == null || !combinedPrMatrix.exists())
            throw new IOException("combined.pr_matrix.tsv missing: " + (combinedPrMatrix == null ? "null" : combinedPrMatrix.getAbsolutePath()));

        resultDir.mkdirs();

        cb.onMessage("DIA post-analysis: reading precursors from combined.parquet");
        DiaNNParquetReader diannReader = new DiaNNParquetReader(combinedParquet.getAbsolutePath(), "KR");
        DiaNNPrecursor[] precursorsQ01 = diannReader.getDiaNNPrecursors();
        if (precursorsQ01 == null) precursorsQ01 = new DiaNNPrecursor[0];

        // Filter and score (matches MetaDiaTask.doExportReportInternal()).
        DiaNNPrecursor[] filtered = filterAndSetRazorProteins(precursorsQ01);

        cb.onMessage("DIA post-analysis: reading peptide matrix and mapping razor proteins");
        DiannPrMatrixReader prReader = new DiannPrMatrixReader(combinedPrMatrix);
        FlashLfqQuanPeptide[] peps = prReader.getMetaPeptides(filtered);
        String[] expNames = prReader.getIntensityTitle();
        if (expNames == null) expNames = new String[0];

        // Build intensity + count maps, and export peptide/protein TSVs.
        HashMap<String, double[]> genomeIntensityMap = new HashMap<>();
        HashMap<String, double[]> proteinIntensityMap = new HashMap<>();
        HashMap<String, Double> proteinTotalIntensityMap = new HashMap<>();
        HashMap<String, Integer> psmCountMap = new HashMap<>();

        HashMap<String, HashSet<String>> totalPepCountMap = new HashMap<>();
        HashMap<String, HashSet<String>> razorPepCountMap = new HashMap<>();

        int[] filePsmCount = new int[expNames.length];
        @SuppressWarnings("unchecked")
        HashSet<String>[] fileSequenceSets = new HashSet[expNames.length];
        for (int i = 0; i < fileSequenceSets.length; i++) fileSequenceSets[i] = new HashSet<>();

        ArrayList<FlashLfqQuanPeptide> peptideList = new ArrayList<>();
        if (peps != null) {
            for (FlashLfqQuanPeptide pep : peps) {
                if (pep == null) continue;
                peptideList.add(pep);
                String stripeSeq = pep.getSequence();
                String razorProtein = pep.getRazorProtein();

                if (stripeSeq != null && razorProtein != null) {
                    razorPepCountMap.computeIfAbsent(razorProtein, k -> new HashSet<>()).add(stripeSeq);
                    if (razorProtein.startsWith(magDb.getIdentifier())) {
                        String genome = razorProtein.split("_")[0];
                        razorPepCountMap.computeIfAbsent(genome, k -> new HashSet<>()).add(stripeSeq);
                    }
                }
            }
        }

        // Score maps for proteins/genomes (S,p_emp,q) — delegated to shared DiaScoring.
        HashMap<String, double[]> genomeDiaStats = DiaScoring.computeGenomeScoreMap(filtered, true);
        HashMap<String, double[]> proteinDiaStats = DiaScoring.computeProteinScoreMap(filtered);

        // All proteins: each peptide's intensity divided among all its proteins (matches MetaDiaMarkerTask).
        HashMap<String, double[]> proteinIntensityMapAll = new HashMap<>();
        HashMap<String, Double> proteinTotalIntenMapAll = new HashMap<>();
        HashMap<String, HashSet<String>> totalPepCountMapAll = new HashMap<>();
        HashMap<String, Integer> psmCountMapAll = new HashMap<>();
        if (peps != null) {
            for (FlashLfqQuanPeptide pep : peps) {
                if (pep == null) continue;
                String[] prots = pep.getProteins();
                if (prots == null || prots.length == 0) continue;
                String stripeSequence = pep.getSequence();
                double[] fullIntensity = pep.getIntensity();
                if (fullIntensity == null) continue;
                int nProts = prots.length;
                double[] divided = new double[fullIntensity.length];
                for (int j = 0; j < fullIntensity.length; j++) {
                    divided[j] = fullIntensity[j] / (double) nProts;
                }
                int[] ms2Count = pep.getMs2Counts();
                int pepMs2Count = 0;
                if (ms2Count != null) {
                    for (int count : ms2Count) pepMs2Count += count;
                }
                for (String pro : prots) {
                    if (pro == null || pro.isBlank()) continue;
                    double[] intensity = proteinIntensityMapAll.get(pro);
                    if (intensity == null) {
                        intensity = new double[divided.length];
                        System.arraycopy(divided, 0, intensity, 0, divided.length);
                        proteinIntensityMapAll.put(pro, intensity);
                        double total = 0.0;
                        for (double v : divided) total += v;
                        proteinTotalIntenMapAll.put(pro, total);
                        HashSet<String> set = new HashSet<>();
                        if (stripeSequence != null) set.add(stripeSequence);
                        totalPepCountMapAll.put(pro, set);
                        psmCountMapAll.put(pro, pepMs2Count);
                    } else {
                        double total = proteinTotalIntenMapAll.getOrDefault(pro, 0.0);
                        for (int k = 0; k < divided.length; k++) {
                            intensity[k] += divided[k];
                            total += divided[k];
                        }
                        proteinTotalIntenMapAll.put(pro, total);
                        if (stripeSequence != null) totalPepCountMapAll.get(pro).add(stripeSequence);
                        psmCountMapAll.put(pro, psmCountMapAll.get(pro) + pepMs2Count);
                    }
                }
            }
        }

        // Aggregate intensities by protein/genome with shared peptide splitting (matches GUI).
        if (peps != null) {
            for (FlashLfqQuanPeptide pep : peps) {
                if (pep == null) continue;
                String[] proteins = pep.getProteins();
                if (proteins == null) continue;

                int count = 0;
                for (String pro : proteins) {
                    if (pro != null && razorPepCountMap.containsKey(pro)) count++;
                }
                if (count <= 0) continue;

                String modSeq = pep.getModSeq();
                int[] idenType = pep.getIdenType();
                if (idenType != null) {
                    for (int j = 0; j < Math.min(fileSequenceSets.length, idenType.length); j++) {
                        if (idenType[j] == 0 && modSeq != null) fileSequenceSets[j].add(modSeq);
                    }
                }

                String stripeSeq = pep.getSequence();
                double[] pepIntensity = new double[pep.getIntensity().length];
                for (int j = 0; j < pepIntensity.length; j++) {
                    pepIntensity[j] = pep.getIntensity()[j] / (double) count;
                }

                int pepMs2Count = 0;
                int[] ms2Count = pep.getMs2Counts();
                if (ms2Count != null) {
                    for (int j = 0; j < Math.min(filePsmCount.length, ms2Count.length); j++) {
                        pepMs2Count += ms2Count[j];
                        filePsmCount[j] += ms2Count[j];
                    }
                }

                for (String pro : proteins) {
                    if (pro == null) continue;
                    if (!razorPepCountMap.containsKey(pro)) continue;

                    psmCountMap.put(pro, psmCountMap.getOrDefault(pro, 0) + pepMs2Count);

                    double[] proInt = proteinIntensityMap.computeIfAbsent(pro, k -> new double[pepIntensity.length]);
                    double total = proteinTotalIntensityMap.getOrDefault(pro, 0.0);
                    for (int k = 0; k < pepIntensity.length; k++) {
                        proInt[k] += pepIntensity[k];
                        total += pepIntensity[k];
                    }
                    proteinTotalIntensityMap.put(pro, total);

                    if (stripeSeq != null) totalPepCountMap.computeIfAbsent(pro, k -> new HashSet<>()).add(stripeSeq);

                    if (pro.startsWith(magDb.getIdentifier())) {
                        String genome = pro.split("_")[0];
                        if (genomeDiaStats.containsKey(genome)) {
                            double[] gInt = genomeIntensityMap.computeIfAbsent(genome, k -> new double[pepIntensity.length]);
                            for (int k = 0; k < pepIntensity.length; k++) gInt[k] += pepIntensity[k];
                            psmCountMap.put(genome, psmCountMap.getOrDefault(genome, 0) + pepMs2Count);
                            if (stripeSeq != null) totalPepCountMap.computeIfAbsent(genome, k -> new HashSet<>()).add(stripeSeq);
                        }
                    }
                }
            }
        }

        File idTsv = new File(resultDir, "id.tsv");
        cb.onMessage("DIA post-analysis: writing " + idTsv.getName());
        try (PrintWriter idWriter = new PrintWriter(idTsv)) {
            idWriter.println("Raw file\tUnique peptide count\tPSM count");
            for (int i = 0; i < expNames.length; i++) {
                idWriter.println(expNames[i] + "\t" + fileSequenceSets[i].size() + "\t" + filePsmCount[i]);
            }
        }

        // Sort proteins by S desc, then total intensity desc (matches GUI).
        String[] proteins = proteinIntensityMap.keySet().toArray(new String[0]);
        Arrays.sort(proteins, (o1, o2) -> {
            double[] s1 = proteinDiaStats.get(o1);
            double[] s2 = proteinDiaStats.get(o2);
            if (s1 == null && s2 == null) return 0;
            if (s1 == null) return 1;
            if (s2 == null) return -1;
            if (s1[IDX_S] > s2[IDX_S]) return -1;
            if (s1[IDX_S] < s2[IDX_S]) return 1;
            double t1 = proteinTotalIntensityMap.getOrDefault(o1, 0.0);
            double t2 = proteinTotalIntensityMap.getOrDefault(o2, 0.0);
            return Double.compare(t2, t1);
        });
        HashMap<String, Integer> proteinIdMap = new HashMap<>();
        for (int i = 0; i < proteins.length; i++) proteinIdMap.put(proteins[i], i);

        File finalPeptides = new File(resultDir, "final_peptides.tsv");
        File finalProteins = new File(resultDir, "final_proteins.tsv");

        cb.onMessage("DIA post-analysis: writing " + finalPeptides.getName());
        HashMap<String, ArrayList<Integer>> proPepIdMap = new HashMap<>();
        HashMap<String, ArrayList<Integer>> proPepRazorIdMap = new HashMap<>();
        writeFinalPeptides(finalPeptides, peptideList, expNames, proteinIdMap, razorPepCountMap, proPepIdMap, proPepRazorIdMap);

        cb.onMessage("DIA post-analysis: writing " + finalProteins.getName());
        writeFinalProteins(finalProteins, proteins, expNames, proteinIntensityMap, proteinTotalIntensityMap, psmCountMap,
                totalPepCountMap, razorPepCountMap, proteinDiaStats, proPepIdMap, proPepRazorIdMap);

        File allProteins = new File(resultDir, "all_proteins.tsv");
        cb.onMessage("DIA post-analysis: writing " + allProteins.getName());
        writeAllProteins(allProteins, peptideList, expNames, proteinIntensityMapAll, proteinTotalIntenMapAll,
                totalPepCountMapAll, razorPepCountMap, psmCountMapAll, proteinDiaStats);

        File finalSummary = new File(resultDir, "final_summary.tsv");
        cb.onMessage("DIA post-analysis: writing " + finalSummary.getName());
        writeFinalSummary(finalSummary, filtered, expNames);

        MetaProtein[] metaProteins = buildMetaProteins(
                proteins, proteinIntensityMap, proteinTotalIntensityMap, psmCountMap, totalPepCountMap, proteinDiaStats);

        // Taxonomy TSVs: compute unique + I/L maps and selection labels from ranked.tsv.
        HashMap<String, HashSet<String>> uniquePepCountMap = new HashMap<>();
        HashMap<String, HashSet<String>> totalPepWithILMap = new HashMap<>();
        HashMap<String, HashSet<String>> razorPepWithILMap = new HashMap<>();
        HashMap<String, HashSet<String>> uniquePepWithILMap = new HashMap<>();
        DiaScoring.fillDiaGenomeUniqueAndILMaps(peps, magDb.getIdentifier(), razorPepCountMap,
                uniquePepCountMap, totalPepWithILMap, razorPepWithILMap, uniquePepWithILMap,
                genomeIntensityMap.keySet());

        File combinedDir = combinedParquet.getParentFile();
        ensureCombinedRankedTsv(combinedDir, filtered, magDb, markerMode, cb);
        File rankedForLabels = resolveRankedTsvForExport(resultDir, combinedDir, markerGenomeRankedTsv);

        HashMap<String, String> genomeSelectionLabel = new HashMap<>();
        if (rankedForLabels != null && rankedForLabels.exists()) {
            try {
                RankedTsvLabelParser.RankedGenomeSelectionResult labels =
                        RankedTsvLabelParser.buildExportSelectionFromRankedTsv(rankedForLabels, genomeIntensityMap.keySet());
                genomeSelectionLabel.putAll(labels.genomeSelectionLabel);
            } catch (Exception ignored) {
            }
        }

        HashMap<String, Double> genomeScoreForTaxonomy = new HashMap<>();
        for (String g : genomeIntensityMap.keySet()) {
            double[] st = genomeDiaStats.get(g);
            double p = st != null && st.length > IDX_P_EMP ? Math.max(1e-300, st[IDX_P_EMP]) : 1.0;
            genomeScoreForTaxonomy.put(g, p);
        }

        File taxonomyDir = new File(resultDir, "taxonomy_analysis");
        File genomeTsv = new File(taxonomyDir, "Genome.tsv");
        File taxaReport = new File(taxonomyDir, "Taxa_report.tsv");

        cb.onMessage("DIA post-analysis: exporting taxonomy TSVs");
        MagTaxonomyTsvExporter.exportGenomeAndTaxaReports(
                taxonomyDir,
                magDb.getOriginalTaxaFile(),
                optionalGtdbNcbiMappingTsv,
                genomeIntensityMap,
                genomeScoreForTaxonomy,
                genomeDiaStats,
                psmCountMap,
                genomeSelectionLabel,
                totalPepCountMap,
                razorPepCountMap,
                uniquePepCountMap,
                totalPepWithILMap,
                razorPepWithILMap,
                uniquePepWithILMap,
                expNames,
                cb
        );

        if (rankedForLabels != null && rankedForLabels.exists()) {
            File rankedCopy = new File(taxonomyDir, "genome_ranked.tsv");
            Files.copy(rankedForLabels.toPath(), rankedCopy.toPath(), StandardCopyOption.REPLACE_EXISTING);
        }

        if (funcDbFile != null && funcDbFile.isFile()) {
            DiaMagFunctionalAnnotationExporter.export(
                    resultDir,
                    magDb,
                    funcDbFile,
                    microDbCatalogueId,
                    quantitativeMode,
                    metaProteins,
                    expNames,
                    cb);
        } else {
            cb.onMessage("DIA functional annotation: skipped (func_def.db not configured)");
        }

        return new Result(idTsv, finalPeptides, finalProteins, taxonomyDir, genomeTsv, taxaReport, genomeIntensityMap, genomeDiaStats);
    }

    private static void ensureCombinedRankedTsv(File combinedDir,
                                                DiaNNPrecursor[] precursors,
                                                MagDbItem magDb,
                                                String markerMode,
                                                ProgressCallback cb) throws Exception {
        if (combinedDir == null) return;
        File rankedTsv = new File(combinedDir, "ranked.tsv");
        if (rankedTsv.exists()) return;
        if (precursors == null || precursors.length == 0) return;
        if (magDb == null || magDb.getCurrentFile() == null) return;

        File genomesTaxaTsv = new File(magDb.getCurrentFile(), "genomes-taxa.tsv");
        if (!genomesTaxaTsv.exists()) return;

        String kneeMode = MarkerGenomeSelector.KNEE_MODE_FOCUS;
        if ("COMPREHENSIVE".equalsIgnoreCase(markerMode)
                || "DEEP".equalsIgnoreCase(markerMode)
                || "DEEP_PER_FILE".equalsIgnoreCase(markerMode)) {
            kneeMode = MarkerGenomeSelector.KNEE_MODE_COMPREHENSIVE;
        }

        if (!combinedDir.exists()) combinedDir.mkdirs();
        File evidenceFile = new File(combinedDir, "evidence.tsv");
        File selectedFile = new File(combinedDir, "selected.tsv");
        cb.onMessage("DIA post-analysis: generating combined/ranked.tsv for Selection_label");
        MarkerGenomeSelector selector = new MarkerGenomeSelector();
        selector.setKneeSelectionMode(kneeMode);
        selector.selectMarkerGenomesWithKnee(
                precursors,
                null,
                genomesTaxaTsv,
                evidenceFile.getAbsolutePath(),
                rankedTsv.getAbsolutePath(),
                selectedFile.getAbsolutePath(),
                MarkerGenomeSelector.UNLIMITED_MAXIMUM_SELECTED,
                true);
    }

    private static File resolveRankedTsvForExport(File resultDir, File combinedDir, File markerGenomeRankedTsv) {
        if (markerGenomeRankedTsv != null && markerGenomeRankedTsv.isFile()) {
            return markerGenomeRankedTsv;
        }
        if (combinedDir != null) {
            File combinedRanked = new File(combinedDir, "ranked.tsv");
            if (combinedRanked.isFile()) return combinedRanked;
        }
        if (resultDir != null) {
            File markerRanked = new File(resultDir, "marker/genome_ranked.tsv");
            if (markerRanked.isFile()) return markerRanked;
        }
        return null;
    }

    private static MetaProtein[] buildMetaProteins(String[] proteins,
                                                   HashMap<String, double[]> proteinIntensityMap,
                                                   HashMap<String, Double> proteinTotalIntensityMap,
                                                   HashMap<String, Integer> psmCountMap,
                                                   HashMap<String, HashSet<String>> totalPepCountMap,
                                                   HashMap<String, double[]> proteinDiaStats) {
        if (proteins == null || proteins.length == 0) return new MetaProtein[0];
        ArrayList<MetaProtein> list = new ArrayList<>();
        for (int i = 0; i < proteins.length; i++) {
            String protein = proteins[i];
            double[] intensity = proteinIntensityMap.get(protein);
            if (intensity == null) intensity = new double[0];
            int pepCount = totalPepCountMap.containsKey(protein) ? totalPepCountMap.get(protein).size() : 0;
            int psm = psmCountMap.getOrDefault(protein, 0);
            double[] scores = proteinDiaStats.get(protein);
            double score = scores != null && scores.length > IDX_S ? scores[IDX_S] : 0.0;
            list.add(new MetaProtein(i, 0, protein, pepCount, psm, score, intensity));
        }
        return list.toArray(new MetaProtein[0]);
    }

    private static void writeFinalPeptides(File out,
                                          ArrayList<FlashLfqQuanPeptide> peptideList,
                                          String[] expNames,
                                          HashMap<String, Integer> proteinIdMap,
                                          HashMap<String, HashSet<String>> razorPepCountMap,
                                          HashMap<String, ArrayList<Integer>> proPepIdMap,
                                          HashMap<String, ArrayList<Integer>> proPepRazorIdMap) throws IOException {
        try (PrintWriter writer = new PrintWriter(out)) {
            StringBuilder titlesb = new StringBuilder();
            titlesb.append("Sequence").append("\t");
            titlesb.append("Base sequence").append("\t");
            titlesb.append("Length").append("\t");
            titlesb.append("Missed cleavages").append("\t");
            titlesb.append("Mass").append("\t");
            titlesb.append("Proteins").append("\t");
            titlesb.append("Charges").append("\t");
            titlesb.append("Qvalue").append("\t");
            titlesb.append("PEP").append("\t");
            for (String exp : expNames) titlesb.append("Identification type ").append(exp).append("\t");
            for (String exp : expNames) titlesb.append("MS/MS count ").append(exp).append("\t");
            titlesb.append("Intensity").append("\t");
            for (String exp : expNames) titlesb.append("Intensity ").append(exp).append("\t");
            titlesb.append("Reverse").append("\t");
            titlesb.append("Potential contaminant").append("\t");
            titlesb.append("id").append("\t");
            titlesb.append("Protein group IDs").append("\t");
            titlesb.append("MS/MS Count");
            writer.println(titlesb);

            int id = 0;
            for (FlashLfqQuanPeptide peptide : peptideList) {
                if (peptide == null) continue;
                String fullSeq = peptide.getModSeq();
                String sequence = peptide.getSequence();
                if (fullSeq == null || sequence == null) continue;

                StringBuilder sb = new StringBuilder();
                sb.append(fullSeq).append("\t");
                sb.append(sequence).append("\t");
                sb.append(sequence.length()).append("\t");
                sb.append(peptide.getMissCleave()).append("\t");
                sb.append(peptide.getMass()).append("\t");

                String[] pros = peptide.getProteins();
                int kept = 0;
                if (pros != null) {
                    for (String pro : pros) {
                        if (pro != null && proteinIdMap.containsKey(pro)) {
                            sb.append(pro).append(";");
                            kept++;
                        }
                    }
                }
                if (kept > 0) sb.deleteCharAt(sb.length() - 1);
                sb.append("\t");

                int[] charge = peptide.getCharges();
                if (charge != null) {
                    StringBuilder chargesb = new StringBuilder();
                    for (int j = 0; j < charge.length; j++) {
                        if (charge[j] == 1) chargesb.append(j).append(";");
                    }
                    if (chargesb.length() > 1) sb.append(chargesb, 0, chargesb.length() - 1).append("\t");
                    else sb.append("2").append("\t");
                } else {
                    sb.append("2").append("\t");
                }

                sb.append(peptide.getScore()).append("\t");
                sb.append(peptide.getPEP()).append("\t");

                int[] idenType = peptide.getIdenType();
                if (idenType != null) {
                    for (int j = 0; j < idenType.length; j++) {
                        int idx = Math.max(0, Math.min(idenType[j], MetaPeptide.idenTypeStrings.length - 1));
                        sb.append(MetaPeptide.idenTypeStrings[idx]).append("\t");
                    }
                } else {
                    for (int i = 0; i < expNames.length; i++) sb.append("\t");
                }

                int totalMsCount = 0;
                int[] ms2Count = peptide.getMs2Counts();
                if (ms2Count != null) {
                    for (int j = 0; j < ms2Count.length; j++) {
                        totalMsCount += ms2Count[j];
                        sb.append(ms2Count[j]).append("\t");
                    }
                } else {
                    for (int i = 0; i < expNames.length; i++) sb.append("0\t");
                }

                double[] intensity = peptide.getIntensity();
                double totalIntensity = 0;
                if (intensity != null) {
                    for (double v : intensity) totalIntensity += v;
                }
                sb.append(BigInteger.valueOf((long) Math.floor(totalIntensity))).append("\t");

                if (intensity != null) {
                    for (double v : intensity) sb.append(BigInteger.valueOf((long) Math.floor(v))).append("\t");
                } else {
                    for (int i = 0; i < expNames.length; i++) sb.append("0\t");
                }

                sb.append("").append("\t"); // reverse
                sb.append("").append("\t"); // contaminant
                sb.append(id).append("\t");

                if (pros != null) {
                    for (String pro : pros) {
                        if (pro != null && proteinIdMap.containsKey(pro)) {
                            sb.append(proteinIdMap.get(pro)).append(";");
                            proPepIdMap.computeIfAbsent(pro, k -> new ArrayList<>()).add(id);

                            HashSet<String> rSet = razorPepCountMap.get(pro);
                            int isRazor = (rSet != null && rSet.contains(sequence)) ? 1 : 0;
                            proPepRazorIdMap.computeIfAbsent(pro, k -> new ArrayList<>()).add(isRazor);
                        }
                    }
                    if (sb.charAt(sb.length() - 1) == ';') sb.deleteCharAt(sb.length() - 1);
                }
                sb.append("\t");
                sb.append(totalMsCount);
                writer.println(sb);
                id++;
            }
        }
    }

    private static void writeAllProteins(File out,
                                         ArrayList<FlashLfqQuanPeptide> peptideList,
                                         String[] expNames,
                                         HashMap<String, double[]> proteinIntensityMapAll,
                                         HashMap<String, Double> proteinTotalIntenMapAll,
                                         HashMap<String, HashSet<String>> totalPepCountMapAll,
                                         HashMap<String, HashSet<String>> razorPepCountMap,
                                         HashMap<String, Integer> psmCountMapAll,
                                         HashMap<String, double[]> proteinScoreMap) throws IOException {
        String[] proteinsAll = proteinTotalIntenMapAll.keySet().toArray(new String[0]);
        Arrays.sort(proteinsAll, Comparator.<String>comparingDouble(
                        p -> proteinScoreMap.containsKey(p) ? -proteinScoreMap.get(p)[IDX_S] : 0.0)
                .thenComparing((p1, p2) -> Double.compare(
                        proteinTotalIntenMapAll.getOrDefault(p2, 0.0),
                        proteinTotalIntenMapAll.getOrDefault(p1, 0.0))));

        HashMap<String, ArrayList<Integer>> proPepIdMapAll = new HashMap<>();
        HashMap<String, ArrayList<Integer>> proPepRazorIdMapAll = new HashMap<>();
        for (int id = 0; id < peptideList.size(); id++) {
            FlashLfqQuanPeptide pep = peptideList.get(id);
            if (pep == null) continue;
            String razorPro = pep.getRazorProtein();
            String[] pros = pep.getProteins();
            if (pros == null) continue;
            for (String pro : pros) {
                if (pro == null || pro.isBlank()) continue;
                proPepIdMapAll.computeIfAbsent(pro, k -> new ArrayList<>()).add(id);
                proPepRazorIdMapAll.computeIfAbsent(pro, k -> new ArrayList<>())
                        .add(pro.equals(razorPro) ? 1 : 0);
            }
        }

        try (PrintWriter writer = new PrintWriter(out)) {
            StringBuilder header = new StringBuilder();
            header.append("Protein IDs").append("\t");
            header.append("Peptide counts (all)").append("\t");
            header.append("Peptide counts (razor)").append("\t");
            header.append("Number of PSMs").append("\t");
            header.append("S").append("\t");
            header.append("p_emp").append("\t");
            header.append("q_value").append("\t");
            header.append("Intensity").append("\t");
            for (String exp : expNames) header.append("Intensity ").append(exp).append("\t");
            header.append("Reverse").append("\t");
            header.append("Potential contaminant").append("\t");
            header.append("id").append("\t");
            header.append("Peptide IDs").append("\t");
            header.append("Peptide is razor");
            writer.println(header);

            DecimalFormat df = FormatTool.getDFE4();
            DecimalFormat df4 = FormatTool.getDF4();
            for (int i = 0; i < proteinsAll.length; i++) {
                String pro = proteinsAll[i];
                int pepAll = totalPepCountMapAll.containsKey(pro) ? totalPepCountMapAll.get(pro).size() : 0;
                int pepRazor = razorPepCountMap.containsKey(pro) ? razorPepCountMap.get(pro).size() : 0;
                StringBuilder sb = new StringBuilder();
                sb.append(pro).append("\t").append(pepAll).append("\t").append(pepRazor).append("\t");
                sb.append(psmCountMapAll.getOrDefault(pro, 0)).append("\t");
                double[] sc = proteinScoreMap.get(pro);
                if (sc != null) {
                    sb.append(df.format(sc[IDX_S])).append("\t");
                    sb.append(df4.format(sc[IDX_P_EMP])).append("\t");
                    sb.append(df4.format(sc[IDX_Q_VALUE])).append("\t");
                } else {
                    sb.append("\t\t\t");
                }
                sb.append(df.format(proteinTotalIntenMapAll.getOrDefault(pro, 0.0))).append("\t");
                double[] intensity = proteinIntensityMapAll.get(pro);
                if (intensity != null) {
                    for (double v : intensity) sb.append(df.format(v)).append("\t");
                } else {
                    for (int j = 0; j < expNames.length; j++) sb.append(df.format(0.0)).append("\t");
                }
                sb.append("\t\t").append(i).append("\t");
                ArrayList<Integer> ids = proPepIdMapAll.get(pro);
                ArrayList<Integer> raz = proPepRazorIdMapAll.get(pro);
                if (ids != null) {
                    for (Integer id : ids) sb.append(id).append(";");
                    if (!ids.isEmpty()) sb.setLength(sb.length() - 1);
                }
                sb.append("\t");
                if (raz != null) {
                    for (Integer r : raz) sb.append(r == 1 ? "true;" : "false;");
                    if (!raz.isEmpty()) sb.setLength(sb.length() - 1);
                }
                writer.println(sb);
            }
        }
    }

    private static void writeFinalSummary(File out, DiaNNPrecursor[] precursors, String[] expNames) throws IOException {
        try (PrintWriter countWriter = new PrintWriter(out)) {
            StringBuilder titlesb = new StringBuilder();
            titlesb.append("Raw file").append("\t");
            titlesb.append("Experiment").append("\t");
            titlesb.append("MS/MS").append("\t");
            titlesb.append("MS/MS Identified").append("\t");
            titlesb.append("MS/MS Identified [%]").append("\t");
            titlesb.append("Peptide Sequences Identified");
            countWriter.println(titlesb);

            HashMap<String, Integer> totalms2Map = new HashMap<>();
            HashMap<String, HashSet<String>> totalUnipepMap = new HashMap<>();
            HashSet<String> totalPepSet = new HashSet<>();
            if (precursors != null) {
                for (DiaNNPrecursor pp : precursors) {
                    if (pp == null) continue;
                    String run = pp.getRunString();
                    String sequence = pp.getSeqString();
                    if (run == null) run = "";
                    if (sequence == null) sequence = "";
                    totalms2Map.put(run, totalms2Map.getOrDefault(run, 0) + 1);
                    totalUnipepMap.computeIfAbsent(run, k -> new HashSet<>()).add(sequence);
                    totalPepSet.add(sequence);
                }
            }

            int totalms2count = 0;
            for (String file : totalUnipepMap.keySet()) {
                StringBuilder sb = new StringBuilder();
                sb.append(file).append("\t").append(file).append("\t");
                int idCounts = totalms2Map.getOrDefault(file, 0);
                int pepCount = totalUnipepMap.get(file).size();
                sb.append(idCounts).append("\t");
                sb.append(idCounts).append("\t");
                sb.append("1.0").append("\t");
                sb.append(pepCount);
                countWriter.println(sb);
                totalms2count += idCounts;
            }

            StringBuilder totalSb = new StringBuilder();
            totalSb.append("Total").append("\t").append("\t");
            totalSb.append(totalms2count).append("\t");
            totalSb.append(totalms2count).append("\t");
            totalSb.append("1.0").append("\t");
            totalSb.append(totalPepSet.size());
            countWriter.print(totalSb);
        }
    }

    private static void writeFinalProteins(File out,
                                          String[] proteins,
                                          String[] expNames,
                                          HashMap<String, double[]> proteinIntensityMap,
                                          HashMap<String, Double> proteinTotalIntensityMap,
                                          HashMap<String, Integer> psmCountMap,
                                          HashMap<String, HashSet<String>> totalPepCountMap,
                                          HashMap<String, HashSet<String>> razorPepCountMap,
                                          HashMap<String, double[]> proteinScoreMap,
                                          HashMap<String, ArrayList<Integer>> proPepIdMap,
                                          HashMap<String, ArrayList<Integer>> proPepRazorIdMap) throws IOException {
        try (PrintWriter writer = new PrintWriter(out)) {
            StringBuilder titlesb = new StringBuilder();
            titlesb.append("Protein IDs").append("\t");
            titlesb.append("Peptide counts (all)").append("\t");
            titlesb.append("Peptide counts (razor)").append("\t");
            titlesb.append("Number of PSMs").append("\t");
            titlesb.append("S").append("\t");
            titlesb.append("p_emp").append("\t");
            titlesb.append("q_value").append("\t");
            titlesb.append("Intensity").append("\t");
            for (String exp : expNames) titlesb.append("Intensity ").append(exp).append("\t");
            titlesb.append("Reverse").append("\t");
            titlesb.append("Potential contaminant").append("\t");
            titlesb.append("id").append("\t");
            titlesb.append("Peptide IDs").append("\t");
            titlesb.append("Peptide is razor");
            writer.println(titlesb);

            DecimalFormat df = FormatTool.getDFE4();
            DecimalFormat df4 = FormatTool.getDF4();
            for (int i = 0; i < proteins.length; i++) {
                String pro = proteins[i];
                StringBuilder sb = new StringBuilder();
                sb.append(pro).append("\t");
                sb.append(totalPepCountMap.containsKey(pro) ? totalPepCountMap.get(pro).size() : 0).append("\t");
                sb.append(razorPepCountMap.containsKey(pro) ? razorPepCountMap.get(pro).size() : 0).append("\t");
                sb.append(psmCountMap.getOrDefault(pro, 0)).append("\t");

                double[] sc = proteinScoreMap.get(pro);
                if (sc != null && sc.length >= 3) {
                    sb.append(df.format(sc[IDX_S])).append("\t");
                    sb.append(df4.format(sc[IDX_P_EMP])).append("\t");
                    sb.append(df4.format(sc[IDX_Q_VALUE])).append("\t");
                } else {
                    sb.append("\t\t\t");
                }

                double totalIntensity = proteinTotalIntensityMap.getOrDefault(pro, 0.0);
                sb.append(df.format(totalIntensity)).append("\t");

                double[] intensity = proteinIntensityMap.get(pro);
                if (intensity != null) {
                    for (double v : intensity) sb.append(df.format(v)).append("\t");
                } else {
                    for (int j = 0; j < expNames.length; j++) sb.append(df.format(0.0)).append("\t");
                }

                sb.append("").append("\t"); // reverse
                sb.append("").append("\t"); // contaminant
                sb.append(i).append("\t");

                ArrayList<Integer> pepIds = proPepIdMap.get(pro);
                ArrayList<Integer> razorFlags = proPepRazorIdMap.get(pro);
                if (pepIds != null) {
                    for (Integer pid : pepIds) sb.append(pid).append(";");
                    sb.deleteCharAt(sb.length() - 1);
                }
                sb.append("\t");
                if (razorFlags != null) {
                    for (Integer f : razorFlags) sb.append(f).append(";");
                    sb.deleteCharAt(sb.length() - 1);
                }
                writer.println(sb);
            }
        }
    }

    private static DiaNNPrecursor[] filterAndSetRazorProteins(DiaNNPrecursor[] precursors) {
        if (precursors == null) return new DiaNNPrecursor[0];
        DiaNNPrecursor[] q01 = precursors;

        // Filter 1: compute genome p_emp to identify genomes with p_emp==1 for protein pruning.
        HashMap<String, double[]> genomeScoreMapForFilter = DiaScoring.computeGenomeScoreMap(q01, true);
        HashSet<String> genomesPempOne = new HashSet<>();
        for (String g : genomeScoreMapForFilter.keySet()) {
            if (genomeScoreMapForFilter.get(g)[IDX_P_EMP] >= 1.0 - 1e-9) genomesPempOne.add(g);
        }

        ArrayList<DiaNNPrecursor> retained = new ArrayList<>();
        for (DiaNNPrecursor pp : q01) {
            if (pp == null) continue;
            if (pp.getGlobalQvalue() >= 0.01) continue;
            String[] proGroups = pp.getProGroups();
            if (proGroups == null || proGroups.length == 0) continue;
            ArrayList<String> kept = new ArrayList<>();
            for (String pro : proGroups) {
                if (pro == null || pro.isEmpty()) continue;
                if (pro.startsWith(FastaManager.shuffle)) { kept.add(pro); continue; }
                String[] parts = pro.split("_");
                String genome = parts.length > 0 ? parts[0] : null;
                if (genome != null && !genomesPempOne.contains(genome)) kept.add(pro);
            }
            if (!kept.isEmpty()) {
                pp.setProGroups(kept.toArray(new String[0]));
                retained.add(pp);
            }
        }

        DiaNNPrecursor[] filtered2 = retained.toArray(new DiaNNPrecursor[0]);
        HashMap<String, double[]> genomeScoreMap = DiaScoring.computeGenomeScoreMap(filtered2, true);
        HashMap<String, double[]> proteinScoreMap = DiaScoring.computeProteinScoreMap(filtered2);

        ArrayList<DiaNNPrecursor> finalList = new ArrayList<>();
        for (DiaNNPrecursor pp : filtered2) {
            if (pp == null) continue;
            String[] pros = pp.getProGroups();
            if (pros == null || pros.length == 0) continue;
            ArrayList<String> kept = new ArrayList<>();
            String razor = "";
            double bestS = -Double.MAX_VALUE;
            for (String pro : pros) {
                double[] sc = proteinScoreMap.get(pro);
                if (sc == null) continue;
                kept.add(pro);
                double S = sc[IDX_S];
                if (S > bestS) { bestS = S; razor = pro; }
            }
            if (!kept.isEmpty()) {
                pp.setProGroups(kept.toArray(new String[0]));
                pp.setRazorProtein(razor);
                finalList.add(pp);
            }
        }

        if (genomeScoreMap.isEmpty()) {
            // keep going; downstream will simply export empty Genome.tsv.
        }

        return finalList.toArray(new DiaNNPrecursor[0]);
    }
}
