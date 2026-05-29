package uk.ac.quadram.metapilot.core.runner;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Set;

/**
 * Parser for final MAG ranked TSV export-time labels (core/mid/tail).
 */
public final class RankedTsvLabelParser {

    private RankedTsvLabelParser() {
    }

    public static final class RankedGenomeSelectionResult {
        public final HashMap<String, String> genomeSelectionLabel = new HashMap<>();
        public final List<String> orderedGenomes = new ArrayList<>();
        public int focusKnee = -1;
        public int comprehensiveKnee = -1;
        public int totalRanked = 0;
        public int privilegeCount = 0;
        public int greedyCount = 0;
    }

    public static RankedGenomeSelectionResult buildExportSelectionFromRankedTsv(File rankedFile, Set<String> allowedGenomes) {
        RankedGenomeSelectionResult result = new RankedGenomeSelectionResult();
        if (rankedFile == null || !rankedFile.isFile()) return result;

        try (BufferedReader reader = new BufferedReader(new FileReader(rankedFile))) {
            String header = reader.readLine();
            if (header == null) return result;
            String[] headerCols = header.split("\t");
            int rankIdx = -1, genomeIdx = -1, greedyGainIdx = -1, reasonIdx = -1;
            int kneeStepIdx = -1, kneeKneedleIdx = -1, kneeCurvIdx = -1, kneePieceIdx = -1, kneeSlopeIdx = -1;
            for (int i = 0; i < headerCols.length; i++) {
                String col = headerCols[i].trim();
                if (col.equals("Rank")) rankIdx = i;
                else if (col.equals("Genome")) genomeIdx = i;
                else if (col.equals("Reason")) reasonIdx = i;
                else if (col.equals("GreedyGain")) greedyGainIdx = i;
                else if (col.equals("KneeStep")) kneeStepIdx = i;
                else if (col.equals("KneeKneedle")) kneeKneedleIdx = i;
                else if (col.equals("KneeCurvature")) kneeCurvIdx = i;
                else if (col.equals("KneePiecewise")) kneePieceIdx = i;
                else if (col.equals("KneeSlope")) kneeSlopeIdx = i;
            }
            int maxColIdx = Math.max(rankIdx, Math.max(genomeIdx,
                    Math.max(greedyGainIdx, Math.max(kneeStepIdx, Math.max(kneeKneedleIdx,
                            Math.max(kneeCurvIdx, Math.max(kneePieceIdx, kneeSlopeIdx)))))));
            if (rankIdx < 0 || genomeIdx < 0 || greedyGainIdx < 0 || kneeStepIdx < 0 || kneeKneedleIdx < 0
                    || kneeCurvIdx < 0 || kneePieceIdx < 0 || kneeSlopeIdx < 0) {
                return result;
            }

            String line;
            boolean kneeValuesRead = false;
            boolean greedyStarted = false;
            int greedyContributionCount = 0;
            while ((line = reader.readLine()) != null) {
                line = line.trim();
                if (line.isEmpty() || line.startsWith("#")) continue;
                String[] parts = line.split("\t");
                if (parts.length <= maxColIdx) continue;
                String genome = parts[genomeIdx].trim();
                if (genome.isEmpty()) continue;
                result.totalRanked++;
                boolean privilege = false;
                if (reasonIdx >= 0 && reasonIdx < parts.length) {
                    String reason = parts[reasonIdx].trim();
                    privilege = reason.startsWith("PRIVILEGE");
                    if (privilege) result.privilegeCount++;
                    else if ("GREEDY".equals(reason)) result.greedyCount++;
                }
                if (allowedGenomes != null && !allowedGenomes.contains(genome)) continue;
                try {
                    Integer.parseInt(parts[rankIdx].trim());
                    if (privilege) {
                        result.genomeSelectionLabel.put(genome, "core");
                        result.orderedGenomes.add(genome);
                        continue;
                    }
                    if (!kneeValuesRead) {
                        int kneeStep = parseIntSafe(parts[kneeStepIdx].trim());
                        int kneeKneedle = parseIntSafe(parts[kneeKneedleIdx].trim());
                        int kneeCurv = parseIntSafe(parts[kneeCurvIdx].trim());
                        int kneePiece = parseIntSafe(parts[kneePieceIdx].trim());
                        int kneeSlope = parseIntSafe(parts[kneeSlopeIdx].trim());
                        int[] focusVals = { kneeKneedle, kneeCurv, kneeSlope };
                        Arrays.sort(focusVals);
                        result.focusKnee = focusVals[1];
                        result.comprehensiveKnee = Math.max(kneeStep, kneePiece);
                        kneeValuesRead = true;
                    }
                    if (result.focusKnee <= 0 || result.comprehensiveKnee <= 0) continue;
                    double greedyGain = parseDoubleSafe(parts[greedyGainIdx].trim());
                    if (!greedyStarted) {
                        if (greedyGain == 0.0) {
                            result.genomeSelectionLabel.put(genome, "core");
                            result.orderedGenomes.add(genome);
                            continue;
                        }
                        greedyStarted = true;
                        greedyContributionCount = 1;
                    } else {
                        if (greedyGain == 0.0) {
                            result.genomeSelectionLabel.put(genome, "tail");
                            result.orderedGenomes.add(genome);
                            continue;
                        }
                        greedyContributionCount++;
                    }
                    if (greedyContributionCount <= result.focusKnee) result.genomeSelectionLabel.put(genome, "core");
                    else if (greedyContributionCount <= result.comprehensiveKnee) result.genomeSelectionLabel.put(genome, "mid");
                    else result.genomeSelectionLabel.put(genome, "tail");
                    result.orderedGenomes.add(genome);
                } catch (NumberFormatException ignored) {
                }
            }

            if (result.greedyCount == 0 && result.totalRanked > 0 && result.orderedGenomes.isEmpty()) {
                try (BufferedReader r2 = new BufferedReader(new FileReader(rankedFile))) {
                    r2.readLine();
                    String l2;
                    while ((l2 = r2.readLine()) != null) {
                        l2 = l2.trim();
                        if (l2.isEmpty() || l2.startsWith("#")) continue;
                        String[] p2 = l2.split("\t");
                        if (p2.length <= genomeIdx) continue;
                        String g = p2[genomeIdx].trim();
                        if (g.isEmpty()) continue;
                        if (allowedGenomes != null && !allowedGenomes.contains(g)) continue;
                        if (!result.genomeSelectionLabel.containsKey(g)) {
                            result.genomeSelectionLabel.put(g, "core");
                            result.orderedGenomes.add(g);
                        }
                    }
                } catch (IOException ignored) {
                }
            }
        } catch (IOException ignored) {
            return result;
        }

        return result;
    }

    private static int parseIntSafe(String s) {
        if (s == null || s.trim().isEmpty()) return 0;
        try { return Integer.parseInt(s.trim()); } catch (NumberFormatException e) { return 0; }
    }

    private static double parseDoubleSafe(String s) {
        if (s == null || s.trim().isEmpty()) return 0.0;
        try { return Double.parseDouble(s.trim()); } catch (NumberFormatException e) { return 0.0; }
    }
}

