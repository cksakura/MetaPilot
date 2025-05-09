package org.metapilot.data;

import org.junit.Before;
import org.junit.Test;
import org.nd4j.linalg.api.ndarray.INDArray;
import static org.junit.Assert.*;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;

public class DataProcessorTest {
    private DataProcessor dataProcessor;
    private Path tempFile;

    @Before
    public void setUp() throws IOException {
        dataProcessor = new DataProcessor();
        
        // Create a temporary test file
        tempFile = Files.createTempFile("test_data", ".tsv");
        try (FileWriter writer = new FileWriter(tempFile.toFile())) {
            writer.write("peptide\tfeature1\tfeature2\tintensity\n");
            writer.write("PEPTIDE1\t1.0\t2.0\t10.0\n");
            writer.write("PEPTIDE2\t3.0\t4.0\t20.0\n");
        }
    }

    @Test
    public void testLoadData() throws IOException {
        dataProcessor.loadData(tempFile.toString());
        
        INDArray features = dataProcessor.getFeatures();
        INDArray targets = dataProcessor.getTargets();
        
        assertNotNull("Features should not be null", features);
        assertNotNull("Targets should not be null", targets);
        assertEquals("Should have 2 samples", 2, features.rows());
        assertEquals("Should have 2 features", 2, features.columns());
        assertEquals("Should have 2 targets", 2, targets.length());
    }

    @Test
    public void testSavePredictions() throws IOException {
        dataProcessor.loadData(tempFile.toString());
        
        double[] predictions = new double[]{15.0, 25.0};
        Path outputFile = Files.createTempFile("predictions", ".tsv");
        dataProcessor.savePredictions(predictions, outputFile.toString());
        
        assertTrue("Output file should exist", Files.exists(outputFile));
        assertTrue("Output file should not be empty", Files.size(outputFile) > 0);
        
        // Clean up
        Files.delete(outputFile);
    }

    @Test
    public void testGetPeptides() throws IOException {
        dataProcessor.loadData(tempFile.toString());
        
        assertEquals("Should have 2 peptides", 2, dataProcessor.getPeptides().size());
        assertEquals("First peptide should match", "PEPTIDE1", dataProcessor.getPeptides().get(0));
        assertEquals("Second peptide should match", "PEPTIDE2", dataProcessor.getPeptides().get(1));
    }
} 