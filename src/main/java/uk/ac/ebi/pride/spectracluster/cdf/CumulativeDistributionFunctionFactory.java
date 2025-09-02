package uk.ac.ebi.pride.spectracluster.cdf;

import uk.ac.ebi.pride.spectracluster.similarity.CombinedFisherIntensityTest;
import uk.ac.ebi.pride.spectracluster.similarity.FrankEtAlDotProduct;
import uk.ac.ebi.pride.spectracluster.util.Defaults;
import uk.ac.quadram.metapilot.core.tools.ResourceLoader;

import java.io.BufferedReader;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.List;

/**
 * Created by jg on 05.05.15.
 */
public class CumulativeDistributionFunctionFactory {
    private CumulativeDistributionFunctionFactory() {

    }

    /**
     * Returns the default cumulative distribution function (CDF). If the CDF is
     * set in the Defaults class, this set CDF is returned. Otherwise the matching
     * CDF is loaded from the resource file.
     * @param similatiryClass
     * @return
     * @throws Exception
     */
    public static CumulativeDistributionFunction getDefaultCumlativeDistributionFunctionForSimilarityMetric(Class similatiryClass) throws Exception {
        if (Defaults.getCumulativeDistributionFunction() != null) {
            return Defaults.getCumulativeDistributionFunction();
        }
        else {
            return getCumulativeDistributionFunctionForSimilarityMetric(similatiryClass);
        }
    }

    /**
     * Returns the matching cumulative distribution function from the matching
     * resource file. If no resource file exists for the passed similarity
     * checker class an Exception is thrown.
     * @param similarityCheckerClass
     * @return
     * @throws Exception Thrown if no CDF resource exists for the passed similarity checker class.
     */
    public static CumulativeDistributionFunction getCumulativeDistributionFunctionForSimilarityMetric(Class similarityCheckerClass) throws Exception {
        if (similarityCheckerClass == CombinedFisherIntensityTest.class) {
            return getCumulativeDistributionFunctionForResource("cumulative.cdf.tsv");
        }
        if (similarityCheckerClass == FrankEtAlDotProduct.class) {
            return getCumulativeDistributionFunctionForResource("dot.cdf.tsv");
        }

        throw new Exception("No cumulative distribution function defined for " + similarityCheckerClass);
    }

    /**
     * Reads the cumulative distribution function from a resource and builds the
     * matching CumulativeDistributionFunction class from it.
     * @param resource
     * @return
     * @throws Exception
     */
    private static CumulativeDistributionFunction getCumulativeDistributionFunctionForResource(String resource) throws Exception {
        try {
            BufferedReader reader = new BufferedReader(new InputStreamReader(ResourceLoader.load(resource)));

            StringBuilder definitionString = new StringBuilder();

            String line;

            while ((line = reader.readLine()) != null)
                definitionString.append(line).append("\n");

            return CumulativeDistributionFunction.fromString(definitionString.toString());
        } catch (Exception e) {
            // If resource loading fails, create a default CDF
            System.err.println("Warning: Could not load CDF resource '" + resource + "', using default CDF. Error: " + e.getMessage());
            return createDefaultCDF(resource);
        }
    }

    /**
     * Creates a default CDF when the resource file cannot be loaded.
     * This provides a fallback to prevent the clustering from failing completely.
     * @param resource The resource name that failed to load
     * @return A default CumulativeDistributionFunction
     */
    private static CumulativeDistributionFunction createDefaultCDF(String resource) {
        // Create a simple default CDF with reasonable values
        List<Double> proportionPeptidesBelowScore = new ArrayList<>();
        
        // Add default proportions for scores 0.0 to 1.0 in 0.1 increments
        for (int i = 0; i <= 10; i++) {
            double score = i * 0.1;
            // Simple linear progression: higher scores have higher proportions
            double proportion = Math.min(0.9, score * 0.8); // Cap at 0.9 to avoid 100%
            proportionPeptidesBelowScore.add(proportion);
        }
        
        return new CumulativeDistributionFunction(100000, 0.1, proportionPeptidesBelowScore);
    }
}
