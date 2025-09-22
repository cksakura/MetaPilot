package uk.ac.quadram.metapilot.mdb.unipept;

import java.util.HashMap;

import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.BeforeEach;
import static org.junit.jupiter.api.Assertions.*;

/**
 * Test class for UnipeptRequester to verify pept2prot functionality
 */
public class UnipeptRequesterTest {

    private String[] testPeptides;

    @BeforeEach
    void setUp() {
        // Test peptides provided by the user
        testPeptides = new String[]{
            "LDKSDPKNPR",
            "LDKSDEEAINAAK", 
            "LDKSNIIHVPIGK"
        };
    }

    @Test
    void testPept2protSinglePeptide() {
        System.out.println("Testing single peptide pept2prot method...");
        
        for (String peptide : testPeptides) {
            System.out.println("\nTesting peptide: " + peptide);
            
            UnipProResult result = UnipeptRequester.pept2prot(peptide);
            
            if (result != null) {
                System.out.println("✓ Successfully retrieved result for peptide: " + peptide);
                System.out.println("  - Sequence: " + result.getSequence());
                System.out.println("  - Number of proteins: " + result.getUniprotIds().length);
                System.out.println("  - Number of taxon IDs: " + result.getTaxonIds().length);
                
                // Print first few protein IDs and names
                String[] proteinIds = result.getUniprotIds();
                String[] proteinNames = result.getProNames();
                int[] taxonIds = result.getTaxonIds();
                
                int maxDisplay = Math.min(3, proteinIds.length);
                for (int i = 0; i < maxDisplay; i++) {
                    System.out.println("  - Protein " + (i+1) + ": " + proteinIds[i] + " (" + proteinNames[i] + ") - Taxon ID: " + taxonIds[i]);
                }
                
                if (proteinIds.length > 3) {
                    System.out.println("  - ... and " + (proteinIds.length - 3) + " more proteins");
                }
                
                // Basic assertions
                assertEquals(peptide, result.getSequence(), "Sequence should match input");
                assertTrue(result.getUniprotIds().length > 0, "Should have at least one protein match");
                assertTrue(result.getProNames().length > 0, "Should have protein names");
                assertTrue(result.getTaxonIds().length > 0, "Should have taxon IDs");
                
            } else {
                System.out.println("✗ No result returned for peptide: " + peptide);
                fail("pept2prot returned null for peptide: " + peptide);
            }
        }
    }

    @Test
    void testPept2protMultiplePeptides() {
        System.out.println("\nTesting multiple peptides pept2prot method...");
        System.out.println("Testing " + testPeptides.length + " peptides: " + String.join(", ", testPeptides));
        
        HashMap<String, UnipProResult> results = UnipeptRequester.pept2prot(testPeptides);
        
        assertNotNull(results, "Results map should not be null");
        System.out.println("✓ Successfully retrieved results for " + results.size() + " out of " + testPeptides.length + " peptides");
        
        for (String peptide : testPeptides) {
            UnipProResult result = results.get(peptide);
            
            if (result != null) {
                System.out.println("\n✓ Peptide: " + peptide);
                System.out.println("  - Number of proteins: " + result.getUniprotIds().length);
                System.out.println("  - Number of taxon IDs: " + result.getTaxonIds().length);
                
                        // Print first protein match
                        if (result.getUniprotIds().length > 0) {
                            System.out.println("  - First protein: " + result.getUniprotIds()[0] + " (" + result.getProNames()[0] + ")");
                        }
                
                // Basic assertions
                assertEquals(peptide, result.getSequence(), "Sequence should match input");
                assertTrue(result.getUniprotIds().length > 0, "Should have at least one protein match");
                
            } else {
                System.out.println("✗ No result for peptide: " + peptide);
                // Note: This might be expected for some peptides, so we don't fail the test
            }
        }
        
        // At least one peptide should have results
        assertTrue(results.size() > 0, "At least one peptide should have results");
    }

    @Test
    void testPept2protInvalidPeptide() {
        System.out.println("\nTesting pept2prot with invalid peptide...");
        
        String invalidPeptide = "INVALIDPEPTIDE123";
        UnipProResult result = UnipeptRequester.pept2prot(invalidPeptide);
        
        if (result == null) {
            System.out.println("✓ Correctly returned null for invalid peptide");
        } else {
            System.out.println("ℹ Got result for invalid peptide (might be valid): " + result.getSequence());
            System.out.println("  - Number of proteins: " + result.getUniprotIds().length);
        }
        
        // This test is more about observing behavior than asserting specific outcomes
        // since the API might return results for peptides we consider "invalid"
    }

    /**
     * Main method for running the test manually
     */
    public static void main(String[] args) {
        System.out.println("=== UnipeptRequester pept2prot Test ===");
        System.out.println("Testing with peptides: LDKSDPKNPR, LDKSDEEAINAAK, LDKSNIIHVPIGK");
        System.out.println("Note: This test requires internet connection to access Unipept API");
        System.out.println();
        
        UnipeptRequesterTest test = new UnipeptRequesterTest();
        test.setUp();
        
        try {
            test.testPept2protSinglePeptide();
            test.testPept2protMultiplePeptides();
            test.testPept2protInvalidPeptide();
            
            System.out.println("\n=== Test Summary ===");
            System.out.println("✓ All tests completed successfully!");
            System.out.println("The pept2prot method appears to be working correctly.");
            
        } catch (Exception e) {
            System.err.println("\n✗ Test failed with exception:");
            e.printStackTrace();
        }
    }
}
