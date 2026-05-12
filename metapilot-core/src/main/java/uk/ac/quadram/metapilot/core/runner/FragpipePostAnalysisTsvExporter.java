package uk.ac.quadram.metapilot.core.runner;

import uk.ac.quadram.metapilot.core.mag.MagInference;
import uk.ac.quadram.metapilot.core.tools.FormatTool;
import uk.ac.quadram.metapilot.core.workflow.ProgressCallback;
import uk.ac.quadram.metapilot.dbSearch.fragpipe.FragPipePSM;
import uk.ac.quadram.metapilot.dbSearch.fragpipe.FragpipePepReader;
import uk.ac.quadram.metapilot.dbSearch.fragpipe.FragpipePeptide;
import uk.ac.quadram.metapilot.task.io.pep.AbstractMetaPeptideReader;
import uk.ac.quadram.metapilot.task.io.pep.MetaPeptide;
import uk.ac.quadram.metapilot.task.io.pep.MetaPSM;
import uk.ac.quadram.metapilot.task.mag.MagDbItem;

import java.io.File;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;

/**
 * Headless post-analysis and TSV export for FragPipe combined outputs.
 *
 * Current scope (no HTML):
 * - Writes final_peptides.tsv and final_proteins.tsv into the workflow result root.
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

        // Build PSM maps
        HashMap<String, MetaPSM> bestPsmByModSeq = new HashMap<>();
        HashMap<String, int[]> pepPsmCountMap = new HashMap<>();
        HashMap<String, int[]> pepChargeMap = new HashMap<>();
        HashMap<String, Integer> psmCountMap = new HashMap<>();

        for (MetaPSM p : psms) {
            if (p == null) continue;
            String modseq = p.getModSeq();
            if (modseq == null) continue;

            // best PSM per peptide by PEP (lower is better)
            if (!bestPsmByModSeq.containsKey(modseq) || p.getPEP() < bestPsmByModSeq.get(modseq).getPEP()) {
                bestPsmByModSeq.put(modseq, p);
            }

            // charge presence map
            int[] charges = pepChargeMap.computeIfAbsent(modseq, k -> new int[7]);
            int c = p.getCharge();
            if (c >= 0 && c < charges.length) charges[c] = 1;

            // per-protein psm counts
            String[] pros = p.getProteins();
            if (pros != null) {
                for (String pro : pros) {
                    if (pro == null) continue;
                    psmCountMap.put(pro, psmCountMap.getOrDefault(pro, 0) + 1);
                }
            }
        }

        // For peptides that have per-file MS2 counts in FragPipe, use MetaPeptide.getMs2Counts()
        for (FragpipePeptide pep : peps) {
            if (pep == null) continue;
            String modseq = pep.getModSeq();
            if (modseq == null) continue;
            int[] counts = pepPsmCountMap.computeIfAbsent(modseq, k -> new int[fileIdMap.size()]);
            int[] sp = pep.getMs2Counts();
            if (sp != null) {
                for (int i = 0; i < Math.min(sp.length, counts.length); i++) {
                    counts[i] = sp[i];
                }
            }
        }

        // Genome scoring + razor protein selection (EXACT strategy as Windows):
        // use DIA-style scoring when PSMs exist, otherwise legacy.
        HashMap<String, Double> genomeScoreMap;
        HashMap<String, double[]> genomeDiaStats = null;
        if (psms != null && psms.length > 0) {
            MagInference.DiaStyleResult dia = MagInference.scoreGenomesDiaStyle(magDb, peps, psms, AbstractMetaPeptideReader.SHUF);
            genomeScoreMap = dia.genomeScoreMap;
            genomeDiaStats = dia.statsForTaxonomy;
        } else {
            genomeScoreMap = MagInference.scoreGenomesLegacy(magDb, peps, AbstractMetaPeptideReader.SHUF);
        }
        String[] razorProteins = MagInference.selectRazorProteins(magDb, peps, genomeScoreMap, AbstractMetaPeptideReader.SHUF);

        // Build intensity maps + peptide count maps (ported from GUI DDAMagPostAnalysisTask).
        HashMap<String, double[]> genomeIntensityMap = new HashMap<>();
        HashMap<String, double[]> proteinIntensityMap = new HashMap<>();
        HashMap<String, Double> proteinTotalIntenMap = new HashMap<>();

        HashMap<String, HashSet<String>> totalPepCountMap = new HashMap<>();
        HashMap<String, HashSet<String>> razorPepCountMap = new HashMap<>();
        HashMap<String, HashSet<String>> uniquePepCountMap = new HashMap<>();

        HashMap<String, HashSet<String>> totalPepWithILMap = new HashMap<>();
        HashMap<String, HashSet<String>> razorPepWithILMap = new HashMap<>();
        HashMap<String, HashSet<String>> uniquePepWithILMap = new HashMap<>();

        ArrayList<MetaPeptide> peptideList = new ArrayList<>();

        // build razor peptide sets for proteins + genomes (for later counting)
        for (int i = 0; i < peps.length; i++) {
            if (razorProteins[i] == null) continue;
            String stripeSeq = peps[i].getSequence();
            if (stripeSeq == null) continue;
            peptideList.add(peps[i]);
            boolean hasIorL = stripeSeq.contains("I") || stripeSeq.contains("L");

            razorPepCountMap.computeIfAbsent(razorProteins[i], k -> new HashSet<>()).add(stripeSeq);
            if (razorProteins[i].startsWith(magDb.getIdentifier())) {
                String genome = razorProteins[i].split("_")[0];
                razorPepCountMap.computeIfAbsent(genome, k -> new HashSet<>()).add(stripeSeq);
                if (hasIorL) {
                    razorPepWithILMap.computeIfAbsent(genome, k -> new HashSet<>()).add(stripeSeq);
                }
            }
        }

        for (int i = 0; i < peps.length; i++) {
            String[] proteins = peps[i].getProteins();
            if (proteins == null) continue;

            int count = 0;
            for (String p : proteins) {
                if (p != null && razorPepCountMap.containsKey(p)) count++;
            }
            if (count <= 0) continue;

            String stripeSeq = peps[i].getSequence();
            if (stripeSeq == null) continue;
            boolean hasIorL = stripeSeq.contains("I") || stripeSeq.contains("L");

            // unique to one genome check (MAG proteins only, same genome)
            String uniqueGenome = null;
            boolean isUnique = true;
            for (String p : proteins) {
                if (p == null) continue;
                if (p.startsWith(magDb.getIdentifier())) {
                    String genome = p.split("_")[0];
                    if (uniqueGenome == null) uniqueGenome = genome;
                    else if (!uniqueGenome.equals(genome)) { isUnique = false; break; }
                } else {
                    isUnique = false;
                    break;
                }
            }

            double[] pepIntensity = peps[i].getIntensity();
            if (pepIntensity == null) pepIntensity = new double[fileIdMap.size()];
            double[] share = new double[pepIntensity.length];
            for (int j = 0; j < share.length; j++) share[j] = pepIntensity[j] / (double) count;

            for (String pro : proteins) {
                if (pro == null) continue;
                if (!razorPepCountMap.containsKey(pro)) continue;

                if (pro.startsWith(magDb.getIdentifier())) {
                    String genome = pro.split("_")[0];
                    if (!genomeScoreMap.containsKey(genome)) continue;

                    double[] gi = genomeIntensityMap.computeIfAbsent(genome, k -> new double[share.length]);
                    for (int k = 0; k < share.length; k++) gi[k] += share[k];

                    double[] pi = proteinIntensityMap.computeIfAbsent(pro, k -> new double[share.length]);
                    double ptot = proteinTotalIntenMap.getOrDefault(pro, 0.0);
                    for (int k = 0; k < share.length; k++) { pi[k] += share[k]; ptot += share[k]; }
                    proteinTotalIntenMap.put(pro, ptot);

                    totalPepCountMap.computeIfAbsent(genome, k -> new HashSet<>()).add(stripeSeq);
                    totalPepCountMap.computeIfAbsent(pro, k -> new HashSet<>()).add(stripeSeq);
                    if (hasIorL) totalPepWithILMap.computeIfAbsent(genome, k -> new HashSet<>()).add(stripeSeq);

                    if (isUnique && uniqueGenome != null && uniqueGenome.equals(genome)) {
                        uniquePepCountMap.computeIfAbsent(genome, k -> new HashSet<>()).add(stripeSeq);
                        if (hasIorL) uniquePepWithILMap.computeIfAbsent(genome, k -> new HashSet<>()).add(stripeSeq);
                    }
                } else {
                    double[] pi = proteinIntensityMap.computeIfAbsent(pro, k -> new double[share.length]);
                    double ptot = proteinTotalIntenMap.getOrDefault(pro, 0.0);
                    for (int k = 0; k < share.length; k++) { pi[k] += share[k]; ptot += share[k]; }
                    proteinTotalIntenMap.put(pro, ptot);
                    totalPepCountMap.computeIfAbsent(pro, k -> new HashSet<>()).add(stripeSeq);
                }
            }
        }

        // Sort proteins by total intensity
        String[] proteins = proteinTotalIntenMap.keySet().toArray(new String[0]);
        Arrays.sort(proteins, Comparator.comparingDouble((String p) -> proteinTotalIntenMap.getOrDefault(p, 0.0)).reversed());

        // Assign numeric IDs to proteins
        HashMap<String, Integer> proteinIdMap = new HashMap<>();
        for (int i = 0; i < proteins.length; i++) {
            proteinIdMap.put(proteins[i], i);
        }

        // Write peptide TSV (schema mirrors MetaPilotMagTask.generatePeptideReport)
        File pepOut = new File(outDir, "final_peptides.tsv");
        cb.onMessage(taskName + ": writing " + pepOut.getName());
        try (PrintWriter writer = new PrintWriter(pepOut)) {
            StringBuilder header = new StringBuilder();
            header.append("Sequence\tBase sequence\tLength\tMissed cleavages\tMass\tProteins\tCharges\tScore\tPEP\t");
            for (String f : fileList) header.append("Identification type ").append(f).append("\t");
            for (String f : fileList) header.append("MS/MS count ").append(f).append("\t");
            header.append("Intensity\t");
            for (String it : intensityTitles) header.append("Intensity ").append(it).append("\t");
            header.append("Reverse\tPotential contaminant\tid\tProtein group IDs\tMS/MS Count");
            writer.println(header);

            int id = 0;
            for (FragpipePeptide pep : peps) {
                if (pep == null) continue;
                String modSeq = pep.getModSeq();
                MetaPSM bestPsm = bestPsmByModSeq.get(modSeq);
                if (bestPsm == null) continue;

                StringBuilder sb = new StringBuilder();
                sb.append(modSeq).append("\t");
                sb.append(pep.getSequence()).append("\t");
                sb.append(pep.getSequence() != null ? pep.getSequence().length() : 0).append("\t");
                sb.append(bestPsm.getMiss()).append("\t");
                sb.append(bestPsm.getCalMass()).append("\t");

                String[] pros = pep.getProteins();
                if (pros != null) {
                    for (String pro : pros) {
                        if (proteinIdMap.containsKey(pro)) sb.append(pro).append(";");
                    }
                    if (sb.length() > 0 && sb.charAt(sb.length() - 1) == ';') sb.deleteCharAt(sb.length() - 1);
                }
                sb.append("\t");

                int[] charge = pepChargeMap.get(modSeq);
                if (charge != null) {
                    StringBuilder chargesb = new StringBuilder();
                    for (int c = 0; c < charge.length; c++) if (charge[c] == 1) chargesb.append(c).append(";");
                    sb.append(chargesb.length() > 0 ? chargesb.substring(0, chargesb.length() - 1) : "2").append("\t");
                } else {
                    sb.append("2\t");
                }

                sb.append(bestPsm.getScore()).append("\t");
                sb.append(bestPsm.getPEP()).append("\t");

                int totalMsCount = 0;
                int[] psmCounts = pepPsmCountMap.get(modSeq);
                double[] inten = pep.getIntensity();
                if (inten == null) inten = new double[fileIdMap.size()];
                if (psmCounts == null) psmCounts = new int[fileIdMap.size()];

                for (int j = 0; j < fileList.size(); j++) {
                    if (psmCounts[j] > 0) {
                        totalMsCount += psmCounts[j];
                        sb.append("By MS/MS\t");
                    } else if (j < inten.length && inten[j] > 0) {
                        sb.append("MBR\t");
                    } else {
                        sb.append("NotDetected\t");
                    }
                }
                for (int j = 0; j < fileList.size(); j++) {
                    sb.append(psmCounts[j]).append("\t");
                }

                double totalIntensity = 0;
                for (double v : inten) totalIntensity += v;
                sb.append((long) Math.floor(totalIntensity)).append("\t");
                for (double v : inten) sb.append((long) Math.floor(v)).append("\t");

                sb.append("\t\t"); // Reverse, Potential contaminant
                sb.append(id).append("\t");

                if (pros != null) {
                    for (String pro : pros) {
                        if (proteinIdMap.containsKey(pro)) sb.append(proteinIdMap.get(pro)).append(";");
                    }
                    if (sb.length() > 0 && sb.charAt(sb.length() - 1) == ';') sb.deleteCharAt(sb.length() - 1);
                }
                sb.append("\t");
                sb.append(totalMsCount);

                writer.println(sb);
                id++;
            }
        }

        // Write protein TSV (schema mirrors MetaPilotMagTask.generateProteinReport)
        File proOut = new File(outDir, "final_proteins.tsv");
        cb.onMessage(taskName + ": writing " + proOut.getName());
        try (PrintWriter writer = new PrintWriter(proOut)) {
            StringBuilder header = new StringBuilder();
            header.append("Protein IDs\tPeptide counts (all)\tPeptide counts (razor)\tNumber of PSMs\tPEP\tIntensity\t");
            for (String it : intensityTitles) header.append("Intensity ").append(it).append("\t");
            header.append("Reverse\tPotential contaminant\tid\tPeptide IDs\tPeptide is razor");
            writer.println(header);

            for (int i = 0; i < proteins.length; i++) {
                String pro = proteins[i];
                double[] inten = proteinIntensityMap.getOrDefault(pro, new double[fileIdMap.size()]);
                double total = proteinTotalIntenMap.getOrDefault(pro, 0.0);
                int psmCount = psmCountMap.getOrDefault(pro, 0);
                HashSet<String> totalPepSet = totalPepCountMap.getOrDefault(pro, new HashSet<>());
                HashSet<String> razorPepSet = razorPepCountMap.getOrDefault(pro, new HashSet<>());

                StringBuilder sb = new StringBuilder();
                sb.append(pro).append("\t");
                sb.append(totalPepSet.size()).append("\t");
                sb.append(razorPepSet.size()).append("\t");
                sb.append(psmCount).append("\t");
                sb.append("0").append("\t"); // PEP placeholder (not computed here)
                sb.append(FormatTool.getDFE4().format(total)).append("\t");
                for (int j = 0; j < inten.length; j++) sb.append(FormatTool.getDFE4().format(inten[j])).append("\t");
                sb.append("\t\t"); // Reverse, Potential contaminant
                sb.append(i).append("\t");
                sb.append("\t"); // Peptide IDs (optional)
                sb.append(""); // Peptide is razor
                writer.println(sb);
            }
        }

        Result res = new Result();
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
        return res;
    }

    public static final class Result {
        public String[] intensityTitles;
        public HashMap<String, double[]> genomeIntensityMap;
        public HashMap<String, Double> genomeScoreMap;
        /** genome -> double[]{S, p_emp, q_value} when DIA-style scoring used; null otherwise */
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

