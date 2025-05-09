package org.metapilot.model;

import org.junit.Before;
import org.junit.Test;
import org.nd4j.linalg.api.ndarray.INDArray;
import org.nd4j.linalg.factory.Nd4j;
import static org.junit.Assert.*;

public class PeptideIntensityPredictorTest {
    private PeptideIntensityPredictor model;
    private INDArray testFeatures;
    private INDArray testTargets;

    @Before
    public void setUp() {
        model = new PeptideIntensityPredictor(15, 128, 64, 0.3);
        
        // Create test data
        testFeatures = Nd4j.rand(100, 15);
        testTargets = Nd4j.rand(100, 1);
    }

    @Test
    public void testModelInitialization() {
        assertNotNull("Model should be initialized", model);
    }

    @Test
    public void testTraining() {
        model.train(testFeatures, testTargets);
        // No exception means training completed successfully
    }

    @Test
    public void testPrediction() {
        model.train(testFeatures, testTargets);
        double[] predictions = model.predict(testFeatures);
        assertEquals("Number of predictions should match number of samples", 
                     testFeatures.rows(), predictions.length);
    }
} 