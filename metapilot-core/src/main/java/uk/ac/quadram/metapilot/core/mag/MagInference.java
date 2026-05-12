package uk.ac.quadram.metapilot.core.mag;

import org.apache.commons.math3.stat.inference.TTest;
import uk.ac.quadram.metapilot.core.math.MathTool;
import uk.ac.quadram.metapilot.task.io.pep.MetaPeptide;
import uk.ac.quadram.metapilot.task.io.pep.MetaPSM;
import uk.ac.quadram.metapilot.task.mag.MagDbItem;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.Set;
import java.util.HashSet;
import java.util.Comparator;

/**
 * Headless MAG inference utilities (ported from GUI workflow logic).
 *
 * Provides:
 * - genome scoring from peptide evidence (legacy method used in GUI as fallback)
 * - razor protein selection using genome scores
 */
public final class MagInference {

    private MagInference() {
    }

    /**
     * Legacy genome scoring (intensity/PEP sum + decoy-based t-test cutoff).
     *
     * Returns genome -> score where smaller is "more significant" (p-value like),
     * matching the downstream "-log10(p)" usage in exports.
     */
    public static HashMap<String, Double> scoreGenomesLegacy(MagDbItem magDb, MetaPeptide[] quanPeps, String decoySymbol) {
        if (magDb == null) throw new IllegalArgumentException("magDb is null");
        if (quanPeps == null) return new HashMap<>();
        if (decoySymbol == null) decoySymbol = "";

        String revIdentifier = decoySymbol + magDb.getIdentifier();
        HashMap<String, Double> genomeIntensityScoreMap = new HashMap<>();
        double[] peptideIntensityScore = new double[quanPeps.length];

        for (int i = 0; i < quanPeps.length; i++) {
            double[] pepIntensity = quanPeps[i].getIntensity();
            double totalIntensity = pepIntensity != null ? MathTool.getTotal(pepIntensity) : 0.0;
            int[] idenType = quanPeps[i].getIdenType();
            int idenCount = 0;
            if (idenType != null) {
                for (int v : idenType) {
                    if (v == MetaPeptide.IDEN_MSMS) idenCount++;
                }
            }
            if (totalIntensity == 0 || idenCount == 0) {
                continue;
            }

            double[] indiPEP = quanPeps[i].getIndividualPEP();
            for (int j = 0; pepIntensity != null && j < pepIntensity.length; j++) {
                if (idenType != null && j < idenType.length && idenType[j] == MetaPeptide.IDEN_MSMS) {
                    double pep = (indiPEP != null && j < indiPEP.length) ? indiPEP[j] : 0.0;
                    if (pep == 0.0) peptideIntensityScore[i] += pepIntensity[j] * 1.0E7;
                    else peptideIntensityScore[i] += pepIntensity[j] / pep;
                }
            }

            String[] proteins = quanPeps[i].getProteins();
            if (proteins == null) continue;
            for (String pro : proteins) {
                if (pro == null) continue;
                if (pro.startsWith(magDb.getIdentifier())) {
                    String genome = pro.split("_")[0];
                    genomeIntensityScoreMap.put(genome, genomeIntensityScoreMap.getOrDefault(genome, 0.0) + peptideIntensityScore[i]);
                } else if (!revIdentifier.isEmpty() && pro.startsWith(revIdentifier)) {
                    String p2 = magDb.getIdentifier() + pro.substring(revIdentifier.length());
                    String genome = decoySymbol + p2.split("_")[0];
                    genomeIntensityScoreMap.put(genome, genomeIntensityScoreMap.getOrDefault(genome, 0.0) + peptideIntensityScore[i]);
                }
            }
        }

        // Build razor genome intensity map
        HashMap<String, Double> razorGenomeIntenMap = new HashMap<>();
        for (int i = 0; i < quanPeps.length; i++) {
            if (peptideIntensityScore[i] <= 0) continue;
            String[] proteins = quanPeps[i].getProteins();
            if (proteins == null) continue;
            String razorGenome = "";
            double razorGenomeIntensity = 0.0;
            for (String pro : proteins) {
                if (pro == null) continue;
                if (pro.startsWith(magDb.getIdentifier())) {
                    String genome = pro.split("_")[0];
                    double gi = genomeIntensityScoreMap.getOrDefault(genome, 0.0);
                    if (gi > razorGenomeIntensity) {
                        razorGenome = genome;
                        razorGenomeIntensity = gi;
                    }
                } else if (!revIdentifier.isEmpty() && pro.startsWith(revIdentifier)) {
                    String p2 = magDb.getIdentifier() + pro.substring(revIdentifier.length());
                    String genome = decoySymbol + p2.split("_")[0];
                    double gi = genomeIntensityScoreMap.getOrDefault(genome, 0.0);
                    if (gi > razorGenomeIntensity) {
                        razorGenome = genome;
                        razorGenomeIntensity = gi;
                    }
                }
            }
            if (!razorGenome.isEmpty()) {
                razorGenomeIntenMap.put(razorGenome, 0.0);
            }
        }

        for (int i = 0; i < quanPeps.length; i++) {
            if (peptideIntensityScore[i] <= 0) continue;
            String[] proteins = quanPeps[i].getProteins();
            if (proteins == null) continue;
            int count = 0;
            for (String pro : proteins) {
                if (pro == null) continue;
                if (pro.startsWith(magDb.getIdentifier())) {
                    String genome = pro.split("_")[0];
                    if (razorGenomeIntenMap.containsKey(genome)) count++;
                } else if (!revIdentifier.isEmpty() && pro.startsWith(revIdentifier)) {
                    String p2 = magDb.getIdentifier() + pro.substring(revIdentifier.length());
                    String genome = decoySymbol + p2.split("_")[0];
                    if (razorGenomeIntenMap.containsKey(genome)) count++;
                }
            }
            if (count > 0) {
                for (String pro : proteins) {
                    if (pro == null) continue;
                    if (pro.startsWith(magDb.getIdentifier())) {
                        String genome = pro.split("_")[0];
                        if (razorGenomeIntenMap.containsKey(genome)) {
                            razorGenomeIntenMap.put(genome, razorGenomeIntenMap.get(genome) + peptideIntensityScore[i]);
                        }
                    } else if (!revIdentifier.isEmpty() && pro.startsWith(revIdentifier)) {
                        String p2 = magDb.getIdentifier() + pro.substring(revIdentifier.length());
                        String genome = decoySymbol + p2.split("_")[0];
                        if (razorGenomeIntenMap.containsKey(genome)) {
                            razorGenomeIntenMap.put(genome, razorGenomeIntenMap.get(genome) + peptideIntensityScore[i]);
                        }
                    }
                }
            }
        }

        List<Map.Entry<String, Double>> sorted = razorGenomeIntenMap.entrySet().stream()
                .sorted(Map.Entry.comparingByValue(Collections.reverseOrder()))
                .toList();
        String[] razorGenomes = sorted.stream().map(Map.Entry::getKey).toArray(String[]::new);

        ArrayList<Double> revLogIntenList = new ArrayList<>();
        for (String g : razorGenomes) {
            if (!revIdentifier.isEmpty() && g.startsWith(revIdentifier)) {
                double intensity = razorGenomeIntenMap.getOrDefault(g, 0.0);
                if (intensity > 0) revLogIntenList.add(Math.log10(intensity));
            }
        }

        HashMap<String, Double> genomeScoreMap = new HashMap<>();
        if (revLogIntenList.isEmpty()) {
            for (String g : razorGenomes) {
                double intensity = razorGenomeIntenMap.getOrDefault(g, 0.0);
                if (intensity == 0) break;
                genomeScoreMap.put(g, 0.0);
            }
            return genomeScoreMap;
        }

        TTest tTest = new TTest();
        double[] revLogIntensities = new double[0];
        for (int i = 2; i < revLogIntenList.size(); i++) {
            revLogIntensities = new double[i];
            for (int j = 0; j < revLogIntensities.length; j++) revLogIntensities[j] = revLogIntenList.get(j);
            double pValue = tTest.tTest(revLogIntensities[0], revLogIntensities);
            if (pValue < 0.01) break;
        }

        if (revLogIntensities.length == 0) {
            for (String g : razorGenomes) {
                double intensity = razorGenomeIntenMap.getOrDefault(g, 0.0);
                if ((!revIdentifier.isEmpty() && g.startsWith(revIdentifier)) || intensity == 0) break;
                genomeScoreMap.put(g, 0.0);
            }
            return genomeScoreMap;
        }

        double[] used = revLogIntensities;
        if (used.length < 6) {
            double[] extended = new double[6];
            for (int i = 0; i < extended.length; i++) extended[i] = used[i % used.length];
            used = extended;
        }

        for (String g : razorGenomes) {
            double intensity = razorGenomeIntenMap.getOrDefault(g, 0.0);
            if ((!revIdentifier.isEmpty() && g.startsWith(revIdentifier)) || intensity == 0) break;
            double logI = Math.log10(intensity);
            double p = tTest.tTest(logI, used);
            genomeScoreMap.put(g, p);
        }

        return genomeScoreMap;
    }

    /**
     * DIA-style genome scoring exactly matching the Windows implementation:
     * - collapse PSM -> peptide (modseq|charge|isDecoy)
     * - e = PSM score (capped), intensityFactor from log2(1+I) normalized into [0.5,1.0]
     * - unique peptides only (shared peptides skipped)
     * - per-genome S = sum(top-k wei)
     * - decoy resampling -> p_emp
     * - BH -> q_value
     * - returns genomeScoreMap = 1 - q_value (Windows behavior)
     * - returns statsForTaxonomy = genome -> {S, p_emp, q_value}
     */
    public static DiaStyleResult scoreGenomesDiaStyle(MagDbItem magDb,
                                                     MetaPeptide[] quanPeps,
                                                     MetaPSM[] psms,
                                                     String decoySymbol) {
        DiaStyleResult out = new DiaStyleResult();
        out.genomeScoreMap = new HashMap<>();
        out.statsForTaxonomy = new HashMap<>();
        if (magDb == null || quanPeps == null || psms == null) {
            return out;
        }
        if (decoySymbol == null) decoySymbol = "";

        String revIdentifier = decoySymbol + magDb.getIdentifier();
        String magId = magDb.getIdentifier();

        HashMap<String, MetaPeptide> modSeqToPeptide = new HashMap<>();
        for (MetaPeptide p : quanPeps) {
            if (p != null && p.getModSeq() != null) {
                modSeqToPeptide.put(p.getModSeq(), p);
            }
        }

        // 1) collapse PSMs
        HashMap<String, DdaPeptideEvidence> peptideMap = new HashMap<>();
        for (MetaPSM psm : psms) {
            if (psm == null) continue;
            String modseq = psm.getModSeq();
            if (modseq == null || modseq.isEmpty()) continue;
            int charge = psm.getCharge();
            double score = psm.getScore();
            String[] pros = psm.getProteins();
            boolean isDecoy = false;
            if (pros != null) {
                for (String pro : pros) {
                    if (pro != null && pro.startsWith(revIdentifier)) {
                        isDecoy = true;
                        break;
                    }
                }
            }
            String key = modseq + "|" + charge + "|" + isDecoy;
            DdaPeptideEvidence ev = peptideMap.get(key);
            if (ev == null) {
                ev = new DdaPeptideEvidence();
                ev.bestScore = score;
                ev.totalIntensity = 0.0;
                ev.isDecoy = isDecoy;
                ev.proteins = new HashSet<>();
                MetaPeptide pep = modSeqToPeptide.get(modseq);
                if (pep != null && pep.getIntensity() != null) {
                    ev.totalIntensity = MathTool.getTotal(pep.getIntensity());
                }
                peptideMap.put(key, ev);
            } else {
                ev.bestScore = Math.max(ev.bestScore, score);
            }
            if (pros != null) {
                for (String pro : pros) {
                    if (pro != null && !pro.isEmpty()) ev.proteins.add(pro);
                }
            }
        }

        // 2) intensity normalization
        double minLog = Double.MAX_VALUE, maxLog = -Double.MAX_VALUE;
        for (DdaPeptideEvidence ev : peptideMap.values()) {
            double logI = Math.log(1.0 + Math.max(0.0, ev.totalIntensity)) / Math.log(2.0);
            minLog = Math.min(minLog, logI);
            maxLog = Math.max(maxLog, logI);
        }
        double logRange = Math.max(1e-10, maxLog - minLog);

        final double cap = 50.0;
        ArrayList<Double> decoyEPool = new ArrayList<>();
        HashMap<String, DdaGenomePeptideEntry> genomePeptides = new HashMap<>();

        int nTargetPeptides = 0;
        int nDecoyPeptides = 0;
        int nUniquePeptides = 0;
        int nSharedPeptides = 0;
        HashMap<Integer, Integer> mDistribution = new HashMap<>();

        for (Map.Entry<String, DdaPeptideEvidence> ent : peptideMap.entrySet()) {
            DdaPeptideEvidence ev = ent.getValue();
            double eRaw = Math.min(cap, Math.max(0.0, ev.bestScore));
            double e = eRaw; // Windows: DDA_GENOME_SCORE_INTENSITY_ONLY=false
            double logI = Math.log(1.0 + Math.max(0.0, ev.totalIntensity)) / Math.log(2.0);
            double norm = (logI - minLog) / logRange;
            double intensityFactor = 0.5 + 0.5 * norm;
            String modseq = ent.getKey().split("\\|", 2)[0];

            if (ev.isDecoy) {
                nDecoyPeptides++;
                decoyEPool.add(e * intensityFactor);
                continue;
            }

            HashSet<String> genomes = new HashSet<>();
            for (String pro : ev.proteins) {
                if (pro == null) continue;
                if (pro.startsWith(revIdentifier)) continue;
                if (pro.startsWith(magId)) {
                    String[] parts = pro.split("_");
                    if (parts.length > 0 && !parts[0].isEmpty()) genomes.add(parts[0]);
                }
            }
            if (genomes.isEmpty()) continue;
            int m = genomes.size();
            boolean uniquePeptidesOnly = true;
            if (m != 1 && uniquePeptidesOnly) {
                nSharedPeptides++;
                Integer prev = mDistribution.get(m);
                mDistribution.put(m, prev == null ? 1 : (prev + 1));
                continue;
            }

            nTargetPeptides++;
            nUniquePeptides++;
            {
                Integer prev = mDistribution.get(1);
                mDistribution.put(1, prev == null ? 1 : (prev + 1));
            }
            double w = (m == 1) ? 1.0 : (1.0 / m);
            double wei = w * e * intensityFactor;

            for (String g : genomes) {
                DdaWeightedScore wps = new DdaWeightedScore();
                wps.seq = modseq;
                wps.e = e;
                wps.w = w;
                wps.intensityFactor = intensityFactor;
                wps.we = wei;
                DdaGenomePeptideEntry ge = genomePeptides.computeIfAbsent(g, k -> new DdaGenomePeptideEntry());
                ge.weightedScores.add(wps);
            }
        }

        // 3) per-genome S
        int k = 5;
        Random rnd = new Random(42);
        int B = 10000;
        HashMap<String, Double> SMap = new HashMap<>();
        HashMap<String, Integer> genomeNMap = new HashMap<>();
        HashMap<String, Integer> genomeTotalPeptidesMap = new HashMap<>();
        for (String g : genomePeptides.keySet()) {
            List<DdaWeightedScore> scores = genomePeptides.get(g).weightedScores;
            genomeTotalPeptidesMap.put(g, scores.size());
            scores.sort((a, b) -> Double.compare(b.we, a.we));
            int n = Math.min(k, scores.size());
            double S = 0.0;
            for (int i = 0; i < n; i++) S += scores.get(i).we;
            SMap.put(g, S);
            genomeNMap.put(g, n);
        }

        // 4) decoy resampling -> p_emp
        if (decoyEPool.isEmpty()) decoyEPool.add(0.0);
        double[] decoyArr = new double[decoyEPool.size()];
        for (int i = 0; i < decoyArr.length; i++) decoyArr[i] = decoyEPool.get(i);
        int decoyLen = decoyArr.length;
        int maxM = 0;
        for (Integer m : genomeTotalPeptidesMap.values()) maxM = Math.max(maxM, m);
        double[] tmp = new double[Math.max(maxM, 1)];

        HashMap<String, Double> pempMap = new HashMap<>();
        for (String g : SMap.keySet()) {
            double Sg = SMap.get(g);
            int M = genomeTotalPeptidesMap.get(g);
            int ng = genomeNMap.get(g);
            int count = 0;
            for (int b = 0; b < B; b++) {
                for (int i = 0; i < M; i++) tmp[i] = decoyArr[rnd.nextInt(decoyLen)];
                Arrays.sort(tmp, 0, M);
                double Sb = 0.0;
                for (int i = M - 1; i >= M - ng && i >= 0; i--) Sb += tmp[i];
                if (Sb >= Sg) count++;
            }
            double pemp = (1.0 + count) / (1.0 + B);
            pempMap.put(g, pemp);
        }

        // 5) BH -> q_value; result = 1 - q_value
        String[] genomes = pempMap.keySet().toArray(new String[0]);
        Arrays.sort(genomes, Comparator.comparingDouble(pempMap::get));
        int nGenomes = genomes.length;
        for (int i = 0; i < nGenomes; i++) {
            String g = genomes[i];
            double pemp = pempMap.get(g);
            double q = Math.min(1.0, pemp * nGenomes / (i + 1));
            out.genomeScoreMap.put(g, 1.0 - q);
            out.statsForTaxonomy.put(g, new double[]{ SMap.get(g), pemp, q });
        }

        // expose some runinfo
        out.k = k;
        out.B = B;
        out.nDecoyPeptides = nDecoyPeptides;
        out.nTargetPeptides = nTargetPeptides;
        out.nUniquePeptides = nUniquePeptides;
        out.nSharedPeptides = nSharedPeptides;
        out.mDistribution = mDistribution;
        out.minLog = minLog;
        out.maxLog = maxLog;
        out.logRange = logRange;
        out.cap = cap;

        return out;
    }

    /**
     * Razor protein selection (ported from GUI).
     */
    public static String[] selectRazorProteins(MagDbItem magDb, MetaPeptide[] quanPeps, HashMap<String, Double> genomeScoreMap, String decoySymbol) {
        if (magDb == null) throw new IllegalArgumentException("magDb is null");
        if (quanPeps == null) return new String[0];
        if (genomeScoreMap == null) genomeScoreMap = new HashMap<>();

        HashMap<String, Double> proteinIntensityMap = new HashMap<>();
        for (MetaPeptide pep : quanPeps) {
            double[] pepIntensity = pep.getIntensity();
            String[] proteins = pep.getProteins();
            if (proteins == null || pepIntensity == null) continue;

            int count = 0;
            for (String pro : proteins) {
                if (pro == null) continue;
                if (decoySymbol != null && !decoySymbol.isEmpty() && pro.startsWith(decoySymbol)) continue;
                if (pro.startsWith(magDb.getIdentifier())) {
                    String genome = pro.split("_")[0];
                    if (genomeScoreMap.containsKey(genome)) count++;
                } else {
                    count++;
                }
            }
            if (count <= 0) continue;

            double intensity = MathTool.getTotal(pepIntensity) / (double) count;
            for (String pro : proteins) {
                if (pro == null) continue;
                if (decoySymbol != null && !decoySymbol.isEmpty() && pro.startsWith(decoySymbol)) continue;
                if (pro.startsWith(magDb.getIdentifier())) {
                    String genome = pro.split("_")[0];
                    if (genomeScoreMap.containsKey(genome)) {
                        proteinIntensityMap.put(pro, proteinIntensityMap.getOrDefault(pro, 0.0) + intensity);
                    }
                } else {
                    proteinIntensityMap.put(pro, proteinIntensityMap.getOrDefault(pro, 0.0) + intensity);
                }
            }
        }

        String[] razorProteins = new String[quanPeps.length];
        for (int i = 0; i < quanPeps.length; i++) {
            String[] proteins = quanPeps[i].getProteins();
            if (proteins == null) continue;

            double bestGenomeScore = Double.MAX_VALUE;
            double bestProteinIntensity = 0.0;
            for (String pro : proteins) {
                if (pro == null) continue;
                if (decoySymbol != null && !decoySymbol.isEmpty() && pro.startsWith(decoySymbol)) continue;

                if (pro.startsWith(magDb.getIdentifier())) {
                    String genome = pro.split("_")[0];
                    if (genomeScoreMap.containsKey(genome)) {
                        double gs = genomeScoreMap.get(genome);
                        double pi = proteinIntensityMap.getOrDefault(pro, 0.0);
                        if (gs < bestGenomeScore || (gs == bestGenomeScore && pi > bestProteinIntensity)) {
                            bestGenomeScore = gs;
                            bestProteinIntensity = pi;
                            razorProteins[i] = pro;
                        }
                    }
                } else {
                    double pi = proteinIntensityMap.getOrDefault(pro, 0.0);
                    if (pi > bestProteinIntensity) {
                        bestProteinIntensity = pi;
                        razorProteins[i] = pro;
                    }
                }
            }
        }
        return razorProteins;
    }

    public static final class DiaStyleResult {
        public HashMap<String, Double> genomeScoreMap;
        /** genome -> double[]{S, p_emp, q_value} */
        public HashMap<String, double[]> statsForTaxonomy;
        // runinfo (optional)
        public int k;
        public int B;
        public int nTargetPeptides;
        public int nDecoyPeptides;
        public int nUniquePeptides;
        public int nSharedPeptides;
        public HashMap<Integer, Integer> mDistribution;
        public double minLog, maxLog, logRange, cap;
    }

    private static final class DdaPeptideEvidence {
        double bestScore;
        double totalIntensity;
        boolean isDecoy;
        Set<String> proteins;
    }

    private static final class DdaWeightedScore {
        String seq;
        double e;
        double w;
        double intensityFactor;
        double we;
    }

    private static final class DdaGenomePeptideEntry {
        ArrayList<DdaWeightedScore> weightedScores = new ArrayList<>();
    }
}

