package uk.ac.quadram.metapilot.task.dia;

import org.junit.jupiter.api.Test;

/**
 * JUnit wrapper around {@link DeepMarkerSingleRawBenchmark} so you can run it via Maven
 * with a correct multi-module classpath.
 *
 * Configure via system properties:
 * - Draw=<rawFile>
 * - Dout=<outDir>
 * - DmagdbJson=<MAGDB_1_0.json>
 * - DcatalogueId=<catalogueId>
 * - DcatalogueVersion=<catalogueVersion>
 * - DgenomesPerRound=<int> (optional, default 5)
 */
public class DeepMarkerSingleRawBenchmarkTest {

	@Test
	public void runBenchmarkOnce() throws Exception {
		String raw = System.getProperty("raw");
		String out = System.getProperty("out");
		String magdbJson = System.getProperty("magdbJson");
		String catalogueId = System.getProperty("catalogueId");
		String catalogueVersion = System.getProperty("catalogueVersion");
		int genomesPerRound = Integer.parseInt(System.getProperty("genomesPerRound", "5"));

		if (raw == null || out == null || magdbJson == null || catalogueId == null || catalogueVersion == null) {
			throw new IllegalArgumentException(
					"Missing system properties. Required: -Draw -Dout -DmagdbJson -DcatalogueId -DcatalogueVersion (optional: -DgenomesPerRound)");
		}

		DeepMarkerSingleRawBenchmark.runOnce(raw, out, magdbJson, catalogueId, catalogueVersion, genomesPerRound);
	}
}

