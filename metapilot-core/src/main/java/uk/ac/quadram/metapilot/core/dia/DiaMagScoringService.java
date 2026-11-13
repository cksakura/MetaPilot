package uk.ac.quadram.metapilot.core.dia;

import uk.ac.quadram.metapilot.core.prodb.FastaManager;
import uk.ac.quadram.metapilot.dbSearch.diann.DiaNNPrecursor;
import uk.ac.quadram.metapilot.quant.flashLFQ.FlashLfqQuanPeptide;

import java.util.*;

/**
 * Canonical, headless DIA MAG scoring utilities shared between headless CLI runners
 * and the Windows GUI (MetaDiaTask, DiaMarkerCoreTask, etc.).
 *
 * <p>All scoring helpers use PEP/CScore + precursor intensity weighting,
 * decoy resampling for empirical p-value, and rank-adjusted q-value calculation.
 * Score arrays are always {@code double[]{S, p_emp, q_value}}.
 */
public final class DiaMagScoringService {

    private DiaMagScoringService() {
    }

    /** Index into score arrays: {@code double[]{S, p_emp, q_value}}. */
    public static final int IDX_S = 0;
    public static final int IDX_P_EMP = 1;
    public static final int IDX_Q_VALUE = 2;

    public static final class FilteredScoringResult {
        public final DiaNNPrecursor[] filteredPrecursors;
        public final HashMap<String, double[]> genomeScoreMap;
        public final HashMap<String, double[]> proteinScoreMap;

        private FilteredScoringResult(DiaNNPrecursor[] filteredPrecursors,
                HashMap<String, double[]> genomeScoreMap,
                HashMap<String, double[]> proteinScoreMap) {
            this.filteredPrecursors = filteredPrecursors;
            this.genomeScoreMap = genomeScoreMap;
            this.proteinScoreMap = proteinScoreMap;
        }
    }

    // ---- Private inner helpers ----

    private static final class PeptideEvidence {
        double bestCScore = 0;
        double bestPEP = 1.0;
        double bestPreQuan = 0;
        double bestMs1Area = 0;
        boolean isDecoy;
        HashSet<String> proteins = new HashSet<>();
    }

    private static final class WeightedPeptideScore {
        double we;
    }

    private static final class GenomePeptideEntry {
        ArrayList<WeightedPeptideScore> weightedScores = new ArrayList<>();
    }

    // ---- Scoring API ----

    /**
     * Apply the canonical DIA post-parquet filter, then compute final genome and
     * protein DIA score maps and assign each retained precursor's razor protein.
     *
     * <p>The filtering strategy is:
     * <ol>
     *   <li>Keep candidate associations for the shared final confidence gate and ranking.</li>
     *   <li>Keep only precursors with {@code Global.Q.Value < 0.01} and at least
     *       one remaining protein.</li>
     *   <li>Recompute genome/protein score maps on the retained precursor set.</li>
     *   <li>Set the razor protein to the retained protein with the highest
     *       protein score {@code S}.</li>
     * </ol>
     *
     * @param precursors DIA precursor array after any caller-specific
     *                   peptide-to-protein/genome mapping has been applied
     * @return filtered precursors plus final genome/protein score maps
     */
    public static FilteredScoringResult filterAndScorePrecursors(DiaNNPrecursor[] precursors) {
        if (precursors == null) precursors = new DiaNNPrecursor[0];

        ArrayList<DiaNNPrecursor> retained = new ArrayList<>();
        for (DiaNNPrecursor pp : precursors) {
            if (pp == null) continue;
            if (pp.getGlobalQvalue() >= 0.01) continue;
            String[] proGroups = pp.getProGroups();
            if (proGroups == null || proGroups.length == 0) continue;

            ArrayList<String> kept = new ArrayList<>();
            for (String pro : proGroups) {
                if (pro == null || pro.isEmpty()) continue;
                if (pro.startsWith(FastaManager.shuffle)) {
                    kept.add(pro);
                    continue;
                }
                String[] parts = pro.split("_");
                String genome = parts.length > 0 ? parts[0] : null;
                if (genome != null) kept.add(pro);
            }

            if (!kept.isEmpty()) {
                pp.setProGroups(kept.toArray(new String[0]));
                retained.add(pp);
            }
        }

        DiaNNPrecursor[] retainedPrecursors = retained.toArray(new DiaNNPrecursor[0]);
        HashMap<String, double[]> proteinScoreMap = computeProteinScoreMap(retainedPrecursors);

        ArrayList<DiaNNPrecursor> finalList = new ArrayList<>();
        for (DiaNNPrecursor pp : retainedPrecursors) {
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
                double s = sc[IDX_S];
                if (s > bestS) {
                    bestS = s;
                    razor = pro;
                }
            }
            if (!kept.isEmpty()) {
                pp.setProGroups(kept.toArray(new String[0]));
                pp.setRazorProtein(razor);
                finalList.add(pp);
            }
        }

        DiaNNPrecursor[] finalPrecursors = finalList.toArray(new DiaNNPrecursor[0]);
        HashMap<String, double[]> finalGenomeScoreMap = computeGenomeScoreMap(finalPrecursors, true);
        HashMap<String, double[]> finalProteinScoreMap = computeProteinScoreMap(finalPrecursors);
        return new FilteredScoringResult(finalPrecursors, finalGenomeScoreMap, finalProteinScoreMap);
    }

    /**
     * Compute per-genome scores {@code {S, p_emp, q_value}} from DIA precursors.
     *
     * <p>Uses top-k (k=5) weighted peptide scores, decoy resampling (B=10000) for
     * empirical p-values, and rank-adjusted q-values.
     *
     * @param precursors        DIA precursor array (may be null or contain nulls)
     * @param uniquePeptidesOnly if true, only peptides mapping to exactly one genome contribute
     * @return map from genome identifier to {@code double[]{S, p_emp, q_value}}
     */
    public static HashMap<String, double[]> computeGenomeScoreMap(DiaNNPrecursor[] precursors,
            boolean uniquePeptidesOnly) {
        HashMap<String, double[]> genomeScoreMap = new HashMap<>();
        if (precursors == null || precursors.length == 0) return genomeScoreMap;

        HashMap<String, PeptideEvidence> peptideMap = buildPeptideEvidenceMap(precursors);

        double[] logRange = computeIntensityLogRange(peptideMap);
        double minLog1 = logRange[0], logRange1 = logRange[2];
        double minLog2 = logRange[1], logRange2 = logRange[3];

        double cap = 6.0;
        final double decoyPepMax = 0.99;
        ArrayList<Double> decoyEPool = new ArrayList<>();
        HashMap<String, GenomePeptideEntry> genomePeptides = new HashMap<>();

        for (PeptideEvidence ev : peptideMap.values()) {
            double e = scoreE(ev.bestCScore, ev.bestPEP, cap);
            double intensityFactor = intensityFactor(ev.bestPreQuan, ev.bestMs1Area,
                    minLog1, logRange1, minLog2, logRange2);

            if (ev.isDecoy) {
                if (ev.bestPEP < decoyPepMax) decoyEPool.add(e * intensityFactor);
                continue;
            }

            HashSet<String> genomes = new HashSet<>();
            for (String pro : ev.proteins) {
                if (pro == null) continue;
                if (pro.startsWith(FastaManager.shuffle)) continue;
                String[] parts = pro.split("_");
                if (parts.length > 0 && !parts[0].isEmpty()) genomes.add(parts[0]);
            }
            if (genomes.isEmpty()) continue;
            int m = genomes.size();
            if (m != 1 && uniquePeptidesOnly) continue;
            double w = (m == 1) ? 1.0 : (1.0 / m);
            double wei = w * e * intensityFactor;

            for (String g : genomes) {
                WeightedPeptideScore wps = new WeightedPeptideScore();
                wps.we = wei;
                genomePeptides.computeIfAbsent(g, k -> new GenomePeptideEntry()).weightedScores.add(wps);
            }
        }

        int k = 5;
        Random rnd = new Random(42);
        int B = 10000;
        HashMap<String, Double> SMap = new HashMap<>();
        HashMap<String, Integer> genomeNMap = new HashMap<>();
        HashMap<String, Integer> genomeTotalPeptidesMap = new HashMap<>();
        for (String g : genomePeptides.keySet()) {
            List<WeightedPeptideScore> scores = genomePeptides.get(g).weightedScores;
            genomeTotalPeptidesMap.put(g, scores.size());
            scores.sort((a, b) -> Double.compare(b.we, a.we));
            int n = Math.min(k, scores.size());
            double S = 0;
            for (int i = 0; i < n; i++) S += scores.get(i).we;
            SMap.put(g, S);
            genomeNMap.put(g, n);
        }

        if (decoyEPool.isEmpty()) decoyEPool.add(0.0);
        double[] decoyArr = decoyEPool.stream().mapToDouble(Double::doubleValue).toArray();
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
                double Sb = 0;
                for (int i = M - 1; i >= M - ng && i >= 0; i--) Sb += tmp[i];
                if (Sb >= Sg) count++;
            }
            pempMap.put(g, (1.0 + count) / (1.0 + B));
        }

        HashMap<String, Double> qValueMap = computeRankAdjustedQValues(pempMap);
        for (String genome : pempMap.keySet()) {
            double pemp = pempMap.get(genome);
            genomeScoreMap.put(genome, new double[]{SMap.get(genome), pemp, qValueMap.get(genome)});
        }
        return genomeScoreMap;
    }

    /**
     * Compute per-protein scores {@code {S, p_emp, q_value}} from DIA precursors.
     *
     * <p>S = sum of all (no top-k) intensity-weighted scores per protein.
     * Decoy resampling and rank-adjusted q-values are applied as for genomes.
     *
     * @param precursors DIA precursor array (may be null or contain nulls)
     * @return map from protein identifier to {@code double[]{S, p_emp, q_value}}
     */
    public static HashMap<String, double[]> computeProteinScoreMap(DiaNNPrecursor[] precursors) {
        HashMap<String, double[]> proteinScoreMap = new HashMap<>();
        if (precursors == null || precursors.length == 0) return proteinScoreMap;

        HashMap<String, PeptideEvidence> peptideMap = buildPeptideEvidenceMap(precursors);

        double[] logRange = computeIntensityLogRange(peptideMap);
        double minLog1 = logRange[0], logRange1 = logRange[2];
        double minLog2 = logRange[1], logRange2 = logRange[3];

        double cap = 6.0;
        final double decoyPepMax = 0.99;
        ArrayList<Double> decoyEPool = new ArrayList<>();
        HashMap<String, ArrayList<Double>> proteinPeptides = new HashMap<>();

        for (PeptideEvidence ev : peptideMap.values()) {
            double e = scoreE(ev.bestCScore, ev.bestPEP, cap);
            double intensityFactor = intensityFactor(ev.bestPreQuan, ev.bestMs1Area,
                    minLog1, logRange1, minLog2, logRange2);
            double wei = e * intensityFactor;

            if (ev.isDecoy) {
                if (ev.bestPEP < decoyPepMax) decoyEPool.add(wei);
                for (String pro : ev.proteins) {
                    if (pro == null || !pro.startsWith(FastaManager.shuffle)) continue;
                    proteinPeptides.computeIfAbsent(pro, k -> new ArrayList<>()).add(wei);
                }
                continue;
            }
            for (String pro : ev.proteins) {
                if (pro == null || pro.startsWith(FastaManager.shuffle)) continue;
                proteinPeptides.computeIfAbsent(pro, k -> new ArrayList<>()).add(wei);
            }
        }

        Random rnd = new Random(42);
        int B = 10000;
        HashMap<String, Double> SMap = new HashMap<>();
        HashMap<String, Integer> proteinNMap = new HashMap<>();
        for (String pro : proteinPeptides.keySet()) {
            ArrayList<Double> weiList = proteinPeptides.get(pro);
            double S = 0;
            for (double ww : weiList) S += ww;
            SMap.put(pro, S);
            proteinNMap.put(pro, weiList.size());
        }

        if (decoyEPool.isEmpty()) decoyEPool.add(0.0);
        double[] decoyArr = decoyEPool.stream().mapToDouble(Double::doubleValue).toArray();
        int decoyLen = decoyArr.length;
        int maxM = 0;
        for (Integer m : proteinNMap.values()) maxM = Math.max(maxM, m);
        double[] tmp = new double[Math.max(maxM, 1)];

        HashMap<String, Double> pempMap = new HashMap<>();
        for (String pro : SMap.keySet()) {
            double Sp = SMap.get(pro);
            int M = proteinNMap.get(pro);
            int count = 0;
            for (int b = 0; b < B; b++) {
                double Sb = 0;
                for (int i = 0; i < M; i++) {
                    tmp[i] = decoyArr[rnd.nextInt(decoyLen)];
                    Sb += tmp[i];
                }
                if (Sb >= Sp) count++;
            }
            pempMap.put(pro, (1.0 + count) / (1.0 + B));
        }

        HashMap<String, Double> qValueMap = computeRankAdjustedQValues(pempMap);
        for (String protein : pempMap.keySet()) {
            double pemp = pempMap.get(protein);
            proteinScoreMap.put(protein, new double[]{SMap.get(protein), pemp, qValueMap.get(protein)});
        }
        return proteinScoreMap;
    }

    /**
     * Fill per-genome unique-peptide and I/L-containing peptide sets.
     *
     * <p>Computes the unique, razor, and total peptide sets (with and without I/L) for
     * each genome in {@code genomeKeys}, populating the six provided maps.
     */
    public static void fillDiaGenomeUniqueAndILMaps(FlashLfqQuanPeptide[] peps,
            String magIdentifier,
            HashMap<String, HashSet<String>> razorPepCountMap,
            HashMap<String, HashSet<String>> uniquePepCountMap,
            HashMap<String, HashSet<String>> totalPepWithILMap,
            HashMap<String, HashSet<String>> razorPepWithILMap,
            HashMap<String, HashSet<String>> uniquePepWithILMap,
            Set<String> genomeKeys) {
        if (peps == null || magIdentifier == null) return;

        for (FlashLfqQuanPeptide pep : peps) {
            if (pep == null) continue;
            String stripeSequence = pep.getSequence();
            String[] proteins = pep.getProteins();
            if (proteins == null) continue;

            int count = 0;
            for (String p : proteins) {
                if (p != null && razorPepCountMap.containsKey(p)) count++;
            }
            if (count == 0) continue;

            boolean hasIorL = stripeSequence != null
                    && (stripeSequence.indexOf('I') >= 0 || stripeSequence.indexOf('L') >= 0);
            String uniqueGenome = null;
            boolean isUnique = true;
            for (String p : proteins) {
                if (p == null || p.isEmpty()) continue;
                if (p.startsWith(magIdentifier)) {
                    String genome = p.split("_")[0];
                    if (uniqueGenome == null) uniqueGenome = genome;
                    else if (!uniqueGenome.equals(genome)) { isUnique = false; break; }
                } else { isUnique = false; break; }
            }

            String razorProtein = pep.getRazorProtein();
            for (String p : proteins) {
                if (p == null || !razorPepCountMap.containsKey(p) || !p.startsWith(magIdentifier)) continue;
                String genome = p.split("_")[0];
                if (!genomeKeys.contains(genome)) continue;
                if (hasIorL) totalPepWithILMap.computeIfAbsent(genome, k -> new HashSet<>()).add(stripeSequence);
                if (isUnique && uniqueGenome != null && uniqueGenome.equals(genome)) {
                    uniquePepCountMap.computeIfAbsent(genome, k -> new HashSet<>()).add(stripeSequence);
                    if (hasIorL) uniquePepWithILMap.computeIfAbsent(genome, k -> new HashSet<>()).add(stripeSequence);
                }
            }
            if (razorProtein != null && razorProtein.startsWith(magIdentifier) && hasIorL) {
                String genome = razorProtein.split("_")[0];
                if (genomeKeys.contains(genome)) {
                    razorPepWithILMap.computeIfAbsent(genome, k -> new HashSet<>()).add(stripeSequence);
                }
            }
        }
    }

    // ---- Private helpers ----

    private static HashMap<String, PeptideEvidence> buildPeptideEvidenceMap(DiaNNPrecursor[] precursors) {
        HashMap<String, PeptideEvidence> peptideMap = new HashMap<>();
        for (DiaNNPrecursor pp : precursors) {
            if (pp == null) continue;
            String seq = pp.getSeqString();
            if (seq == null || seq.isEmpty()) continue;
            String key = seq + "|" + pp.getCharge() + "|" + pp.getIsDecoy();
            PeptideEvidence ev = peptideMap.get(key);
            if (ev == null) {
                ev = new PeptideEvidence();
                ev.bestCScore = pp.getCscore();
                ev.bestPEP = pp.getPEP();
                ev.bestPreQuan = pp.getPreQuan();
                ev.bestMs1Area = pp.getMs1Area();
                ev.isDecoy = pp.getIsDecoy();
                peptideMap.put(key, ev);
            } else {
                ev.bestCScore = Math.max(ev.bestCScore, pp.getCscore());
                ev.bestPEP = Math.min(ev.bestPEP, pp.getPEP());
                ev.bestPreQuan = Math.max(ev.bestPreQuan, pp.getPreQuan());
                ev.bestMs1Area = Math.max(ev.bestMs1Area, pp.getMs1Area());
            }
            String[] proGroups = pp.getProGroups();
            if (proGroups != null) {
                for (String pro : proGroups) if (pro != null && !pro.isEmpty()) ev.proteins.add(pro);
            }
        }
        return peptideMap;
    }

    /**
     * Returns {@code [minLog1, minLog2, logRange1, logRange2]} for intensity normalisation.
     */
    private static double[] computeIntensityLogRange(HashMap<String, PeptideEvidence> peptideMap) {
        double minLog1 = Double.MAX_VALUE, maxLog1 = -Double.MAX_VALUE;
        double minLog2 = Double.MAX_VALUE, maxLog2 = -Double.MAX_VALUE;
        for (PeptideEvidence ev : peptideMap.values()) {
            double l1 = Math.log(1.0 + Math.max(0, ev.bestPreQuan)) / Math.log(2.0);
            double l2 = Math.log(1.0 + Math.max(0, ev.bestMs1Area)) / Math.log(2.0);
            minLog1 = Math.min(minLog1, l1);
            maxLog1 = Math.max(maxLog1, l1);
            minLog2 = Math.min(minLog2, l2);
            maxLog2 = Math.max(maxLog2, l2);
        }
        return new double[]{minLog1, minLog2,
                Math.max(1e-10, maxLog1 - minLog1),
                Math.max(1e-10, maxLog2 - minLog2)};
    }

    private static double scoreE(double cscore, double pep, double cap) {
        return cscore > 0
                ? Math.min(cap, Math.max(0, cscore))
                : Math.min(cap, Math.max(0, -Math.log10(Math.max(pep, 1e-10))));
    }

    private static HashMap<String, Double> computeRankAdjustedQValues(HashMap<String, Double> pValueMap) {
        HashMap<String, Double> qValueMap = new HashMap<>();
        if (pValueMap == null || pValueMap.isEmpty()) return qValueMap;

        String[] ids = pValueMap.keySet().toArray(new String[0]);
        Arrays.sort(ids, Comparator.comparingDouble(pValueMap::get));
        int n = ids.length;
        double minAdjusted = 1.0;
        for (int i = n - 1; i >= 0; i--) {
            double p = pValueMap.get(ids[i]);
            double adjusted = Math.min(1.0, p * n / (i + 1));
            minAdjusted = Math.min(minAdjusted, adjusted);
            qValueMap.put(ids[i], minAdjusted);
        }
        return qValueMap;
    }

    private static double intensityFactor(double preQuan, double ms1Area,
            double minLog1, double logRange1, double minLog2, double logRange2) {
        double l1 = Math.log(1.0 + Math.max(0, preQuan)) / Math.log(2.0);
        double l2 = Math.log(1.0 + Math.max(0, ms1Area)) / Math.log(2.0);
        double n1 = (l1 - minLog1) / logRange1;
        double n2 = (l2 - minLog2) / logRange2;
        return 0.5 + 0.25 * n1 + 0.25 * n2;
    }
}
