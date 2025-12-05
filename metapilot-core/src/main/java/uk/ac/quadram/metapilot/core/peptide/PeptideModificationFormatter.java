package uk.ac.quadram.metapilot.core.peptide;

import java.util.ArrayList;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * Utility class for formatting modified peptide sequences
 * 
 * @author Assistant
 */
public class PeptideModificationFormatter {

    /**
     * Formats a peptide sequence with modifications
     * 
     * @param sequence The base peptide sequence (e.g., "MKDLFERSVPMQETAQTK")
     * @param modInfo The modification information (e.g., "N-term(42.0106),11M(15.9949)")
     * @return The formatted modified sequence (e.g., "n[42.0106]MKDLFERSVPM[15.9949]QETAQTK")
     */
    public static String getModifiedPeptideSequence(String sequence, String modInfo) {
        if (sequence == null || sequence.isEmpty()) {
            return "";
        }
        
        if (modInfo == null || modInfo.trim().isEmpty()) {
            return sequence;
        }
        
        // Parse modifications
        List<Modification> modifications = parseModifications(modInfo);
        
        // Build the modified sequence
        StringBuilder result = new StringBuilder();
        
        // Handle N-terminal modification
        String nTermMod = null;
        for (Modification mod : modifications) {
            if (mod.isNTerminal()) {
                nTermMod = mod.getMass();
                break;
            }
        }
        
        if (nTermMod != null) {
            result.append("n[").append(nTermMod).append("]");
        }
        
        // Process each amino acid in the sequence
        for (int i = 0; i < sequence.length(); i++) {
            char aminoAcid = sequence.charAt(i);
            
            // Check if this position has a modification
            String modMass = null;
            for (Modification mod : modifications) {
                if (mod.getPosition() == i && mod.getAminoAcid() == aminoAcid) {
                    modMass = mod.getMass();
                    break;
                }
            }
            
            // Add the amino acid with or without modification
            if (modMass != null) {
                result.append(aminoAcid).append("[").append(modMass).append("]");
            } else {
                result.append(aminoAcid);
            }
        }
        
        return result.toString();
    }
    
    /**
     * Parses modification information string into Modification objects
     * 
     * @param modInfo The modification string (e.g., "N-term(42.0106),11M(15.9949)")
     * @return List of Modification objects
     */
    private static List<Modification> parseModifications(String modInfo) {
        List<Modification> modifications = new ArrayList<>();
        
        // Split by comma to get individual modifications
        String[] modParts = modInfo.split(",");
        
        for (String modPart : modParts) {
            modPart = modPart.trim();
            
            // Handle N-terminal modification
            if (modPart.startsWith("N-term(")) {
                Pattern pattern = Pattern.compile("N-term\\(([0-9.]+)\\)");
                Matcher matcher = pattern.matcher(modPart);
                if (matcher.find()) {
                    modifications.add(new Modification(true, -1, '\0', matcher.group(1)));
                }
            }
            // Handle position-specific modifications (e.g., "11M(15.9949)")
            else {
                Pattern pattern = Pattern.compile("(\\d+)([A-Z])\\(([0-9.]+)\\)");
                Matcher matcher = pattern.matcher(modPart);
                 if (matcher.find()) {
                     int position = Integer.parseInt(matcher.group(1)) - 1; // Convert to 0-based indexing
                     char aminoAcid = matcher.group(2).charAt(0);
                     String mass = matcher.group(3);
                     modifications.add(new Modification(false, position, aminoAcid, mass));
                 }
            }
        }
        
        return modifications;
    }
    
    /**
     * Inner class to represent a modification
     */
    private static class Modification {
        private boolean nTerminal;
        private int position;
        private char aminoAcid;
        private String mass;
        
        public Modification(boolean nTerminal, int position, char aminoAcid, String mass) {
            this.nTerminal = nTerminal;
            this.position = position;
            this.aminoAcid = aminoAcid;
            this.mass = mass;
        }
        
        public boolean isNTerminal() {
            return nTerminal;
        }
        
        public int getPosition() {
            return position;
        }
        
        public char getAminoAcid() {
            return aminoAcid;
        }
        
        public String getMass() {
            return mass;
        }
    }
}
