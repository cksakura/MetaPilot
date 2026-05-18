package uk.ac.quadram.metapilot.core.runner;

import uk.ac.quadram.metapilot.core.mag.DdaMagPostAnalysisService;
import uk.ac.quadram.metapilot.core.workflow.ProgressCallback;
import uk.ac.quadram.metapilot.dbSearch.fragpipe.FragPipePSM;
import uk.ac.quadram.metapilot.dbSearch.fragpipe.FragpipePepReader;
import uk.ac.quadram.metapilot.dbSearch.fragpipe.FragpipePeptide;
import uk.ac.quadram.metapilot.task.io.pep.AbstractMetaPeptideReader;
import uk.ac.quadram.metapilot.task.io.pep.MetaPSM;
import uk.ac.quadram.metapilot.task.io.pro.MetaProtein;
import uk.ac.quadram.metapilot.task.mag.MagDbItem;
import uk.ac.quadram.metapilot.task.par.MetaConstants;

import java.io.File;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;

/**
 * Headless post-analysis and TSV export for FragPipe combined outputs.
 *
 * Current scope:
 * - Writes id.tsv, final_summary.tsv, final_peptides.tsv, and razor_proteins.tsv.
 *
 * This mirrors the GUI post-analysis maps (genome/protein intensity + peptide counts),
 * without any JavaFX dependencies.
 */
public final class FragpipePostAnalysisTsvExporter {

    private FragpipePostAnalysisTsvExporter() {
    }

    public static Result exportFinalTsvs(MagDbItem magDb,
                                         File combinedModifiedPeptideTsv,
                                         File outDir,
                                         String taskName,
                                         ProgressCallback cb) throws Exception {
        return exportFinalTsvs(magDb, combinedModifiedPeptideTsv, outDir, null, taskName, cb);
    }

    public static Result exportFinalTsvs(MagDbItem magDb,
                                         File combinedModifiedPeptideTsv,
                                         File outDir,
                                         File spectraDir,
                                         String taskName,
                                         ProgressCallback cb) throws Exception {
        if (cb == null) cb = new ProgressCallback.NoOp();
        if (magDb == null) throw new IllegalArgumentException("magDb is null");
        if (combinedModifiedPeptideTsv == null || !combinedModifiedPeptideTsv.exists()) {
            throw new IllegalArgumentException("combined_modified_peptide.tsv not found: " + combinedModifiedPeptideTsv);
        }
        if (outDir == null) throw new IllegalArgumentException("outDir is null");
        outDir.mkdirs();

        cb.onMessage(taskName + ": loading combined peptides/PSMs ...");
        FragpipePepReader reader = new FragpipePepReader(combinedModifiedPeptideTsv);
        FragpipePeptide[] peps = reader.getMetaPeptides();
        FragPipePSM[] psms = reader.getPSMs();
        String[] intensityTitles = reader.getIntensityTitle();

        if (peps == null) peps = new FragpipePeptide[0];
        if (psms == null) psms = new FragPipePSM[0];
        if (intensityTitles == null) intensityTitles = new String[0];

        cb.onMessage(taskName + ": peptides=" + peps.length + ", psms=" + psms.length);

        // Build file list from intensity titles (FragPipe column prefixes)
        ArrayList<String> fileList = new ArrayList<>();
        HashMap<String, Integer> fileIdMap = new HashMap<>();
        for (String t : intensityTitles) {
            if (t == null || t.isBlank()) continue;
            String key = t.trim();
            fileIdMap.put(key, fileIdMap.size());
            fileList.add(key);
        }

        DdaMagPostAnalysisService.Result analysis = DdaMagPostAnalysisService.analyze(
                new DdaMagPostAnalysisService.Input(peps, psms, fileIdMap, fileList,
                        AbstractMetaPeptideReader.SHUF, magDb.getIdentifier(),
                        FragpipePostAnalysisTsvExporter::resolveFileIndex, true));
        HashMap<String, MetaPSM> bestPsmByModSeq = analysis.psmMap;
        HashMap<String, int[]> pepPsmCountMap = analysis.pepPsmCountMap;
        HashMap<String, int[]> pepChargeMap = analysis.pepChargeMap;
        HashMap<String, Integer> psmCountMap = analysis.psmCountMap;
        int[] filePsmCount = analysis.filePsmCount;
        HashSet<String>[] fileSequenceSets = analysis.fileSequenceSets;
        HashMap<String, Double> genomeScoreMap = analysis.genomeScoreMap;
        HashMap<String, double[]> genomeDiaStats = analysis.genomeDiaStats;
        HashMap<String, double[]> genomeIntensityMap = analysis.genomeIntensityMap;
        HashMap<String, double[]> proteinIntensityMap = analysis.proteinIntensityMap;
        HashMap<String, Double> proteinTotalIntenMap = analysis.proteinTotalIntenMap;
        HashMap<String, HashSet<String>> totalPepCountMap = analysis.totalPepCountMap;
        HashMap<String, HashSet<String>> razorPepCountMap = analysis.razorPepCountMap;
        HashMap<String, HashSet<String>> uniquePepCountMap = analysis.uniquePepCountMap;
        HashMap<String, HashSet<String>> totalPepWithILMap = analysis.totalPepWithILMap;
        HashMap<String, HashSet<String>> razorPepWithILMap = analysis.razorPepWithILMap;
        HashMap<String, HashSet<String>> uniquePepWithILMap = analysis.uniquePepWithILMap;
        String[] proteins = analysis.proteins;
        HashMap<String, Integer> proteinIdMap = analysis.proteinIdMap;
        HashMap<String, String> razorProteinMap = buildRazorProteinMap(peps, analysis.razorProteins);
        MetaProtein[] metaProteins = analysis.metaProteins;
        cb.onMessage(taskName + ": genome scoring="
                + (analysis.genomeScoring.diaStyle ? "PSM-score/top-k/q-value" : "legacy fallback")
                + ", genomes=" + genomeScoreMap.size());
        File idOut = new File(outDir, "combined/id.tsv");
        File summaryOut = new File(outDir, "final_summary.tsv");
        File pepOut = new File(outDir, "final_peptides.tsv");
        File proOut = new File(outDir, "razor_proteins.tsv");
        File allProOut = new File(outDir, "all_proteins.tsv");
        cb.onMessage(taskName + ": writing " + idOut.getName());
        DdaMagReportExportService.writeIdentificationSummary(new File(outDir, "combined"), outDir, spectraDir,
                fileIdMap, filePsmCount, fileSequenceSets, null);
        cb.onMessage(taskName + ": writing " + pepOut.getName());
        DdaMagReportExportService.writePeptideReport(pepOut, analysis.peptideList, bestPsmByModSeq,
                pepPsmCountMap, pepChargeMap, proteinIdMap, fileList, MetaConstants.labelFree, intensityTitles,
                razorProteinMap);
        cb.onMessage(taskName + ": writing " + proOut.getName());
        DdaMagReportExportService.writeProteinReport(proOut, proteins, proteinIntensityMap, proteinTotalIntenMap,
                totalPepCountMap, razorPepCountMap, psmCountMap, MetaConstants.labelFree, intensityTitles,
                analysis.peptideList, analysis.allProteinScoreMap);
        cb.onMessage(taskName + ": writing " + allProOut.getName());
        DdaMagReportExportService.writeAllProteinReport(allProOut, analysis.allProteins,
                analysis.allProteinIntensityMap, analysis.allProteinTotalIntenMap, analysis.allProteinPepCountMap,
                razorPepCountMap, analysis.allProteinPsmCountMap, MetaConstants.labelFree, intensityTitles,
                analysis.peptideList, analysis.allProteinScoreMap);
        Result res = new Result();
        res.idTsv = idOut;
        res.finalSummaryTsv = summaryOut;
        res.intensityTitles = intensityTitles;
        res.genomeIntensityMap = genomeIntensityMap;
        res.genomeScoreMap = genomeScoreMap;
        res.genomeDiaStats = genomeDiaStats;
        res.psmCountMap = psmCountMap;
        res.totalPepCountMap = totalPepCountMap;
        res.razorPepCountMap = razorPepCountMap;
        res.uniquePepCountMap = uniquePepCountMap;
        res.totalPepWithILMap = totalPepWithILMap;
        res.razorPepWithILMap = razorPepWithILMap;
        res.uniquePepWithILMap = uniquePepWithILMap;
        res.metaProteins = metaProteins;
        return res;
    }

    private static HashMap<String, String> buildRazorProteinMap(FragpipePeptide[] peptides, String[] razorProteins) {
        HashMap<String, String> map = new HashMap<>();
        if (peptides == null || razorProteins == null) {
            return map;
        }
        for (int i = 0; i < peptides.length && i < razorProteins.length; i++) {
            if (peptides[i] == null || razorProteins[i] == null) {
                continue;
            }
            if (peptides[i].getModSeq() != null) {
                map.put(peptides[i].getModSeq(), razorProteins[i]);
            }
            if (peptides[i].getSequence() != null) {
                map.put(peptides[i].getSequence(), razorProteins[i]);
            }
        }
        return map;
    }

    /**
     * @deprecated Identification summary output is now written by
     *             {@link DdaMagReportExportService#writeIdentificationSummary(File, File, File, HashMap, int[], HashSet[], Map)}.
     */
    @Deprecated
    private static void writeIdentificationTsv(File out,
                                               ArrayList<String> fileList,
                                               HashSet<String>[] fileSequenceSets,
                                               int[] filePsmCount) throws Exception {
        if (out.getParentFile() != null) out.getParentFile().mkdirs();
        try (PrintWriter writer = new PrintWriter(out)) {
            writer.println("Raw file\tUnique peptide count\tPSM count");
            for (int i = 0; i < fileList.size(); i++) {
                int unique = i < fileSequenceSets.length && fileSequenceSets[i] != null ? fileSequenceSets[i].size() : 0;
                int psm = i < filePsmCount.length ? filePsmCount[i] : 0;
                writer.println(fileList.get(i) + "\t" + unique + "\t" + psm);
            }
        }
    }

    /**
     * @deprecated Identification summary output is now written by
     *             {@link DdaMagReportExportService#writeIdentificationSummary(File, File, File, HashMap, int[], HashSet[], Map)}.
     */
    @Deprecated
    private static void writeFinalSummary(File out,
                                          ArrayList<String> fileList,
                                          HashSet<String>[] fileSequenceSets,
                                          int[] filePsmCount) throws Exception {
        if (out.getParentFile() != null) out.getParentFile().mkdirs();
        try (PrintWriter writer = new PrintWriter(out)) {
            writer.println("Raw file\tExperiment\tMS/MS\tMS/MS Identified\tMS/MS Identified [%]\tPeptide Sequences Identified");

            int totalPsmCount = 0;
            int totalUniquePeptideCount = 0;
            for (int i = 0; i < fileList.size(); i++) {
                int psm = i < filePsmCount.length ? filePsmCount[i] : 0;
                int unique = i < fileSequenceSets.length && fileSequenceSets[i] != null ? fileSequenceSets[i].size() : 0;
                writer.println(fileList.get(i) + "\t" + fileList.get(i) + "\t"
                        + psm + "\t" + psm + "\t1.0\t" + unique);
                totalPsmCount += psm;
                totalUniquePeptideCount += unique;
            }

            writer.print("Total\t\t" + totalPsmCount + "\t" + totalPsmCount + "\t1.0\t" + totalUniquePeptideCount);
        }
    }

    private static int resolveFileIndex(String fileName, Map<String, Integer> fileIdMap, List<String> fileList) {
        if (fileName == null || fileName.isBlank()) return -1;
        String key = fileName.trim();
        Integer direct = fileIdMap.get(key);
        if (direct != null) return direct;

        String base = baseNameNoExt(key);
        for (int i = 0; i < fileList.size(); i++) {
            String candidate = fileList.get(i);
            if (candidate == null) continue;
            if (key.equals(candidate.trim())) return i;
            String candidateBase = baseNameNoExt(candidate.trim());
            if (base.equals(candidateBase)) return i;
        }
        return -1;
    }

    private static String baseNameNoExt(String path) {
        String name = new File(path).getName();
        int dot = name.lastIndexOf('.');
        return dot > 0 ? name.substring(0, dot) : name;
    }

    public static final class Result {
        public File idTsv;
        public File finalSummaryTsv;
        public String[] intensityTitles;
        public MetaProtein[] metaProteins;
        public HashMap<String, double[]> genomeIntensityMap;
        public HashMap<String, Double> genomeScoreMap;
        /** genome -> double[]{S, p_emp, q_value} when PSM-score/top-k scoring is used; null otherwise */
        public HashMap<String, double[]> genomeDiaStats;
        public HashMap<String, Integer> psmCountMap;
        public HashMap<String, HashSet<String>> totalPepCountMap;
        public HashMap<String, HashSet<String>> razorPepCountMap;
        public HashMap<String, HashSet<String>> uniquePepCountMap;
        public HashMap<String, HashSet<String>> totalPepWithILMap;
        public HashMap<String, HashSet<String>> razorPepWithILMap;
        public HashMap<String, HashSet<String>> uniquePepWithILMap;
    }
}

