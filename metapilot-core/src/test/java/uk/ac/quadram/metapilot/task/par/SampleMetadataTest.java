package uk.ac.quadram.metapilot.task.par;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.nio.file.Files;
import java.nio.file.Path;
import java.util.List;

import org.json.JSONObject;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;

import uk.ac.quadram.metapilot.task.dda.mag.par.DDAMagFragPipeParameterIO;
import uk.ac.quadram.metapilot.task.dda.mag.par.DDAMagFragPipeParameters;

class SampleMetadataTest {

	@TempDir
	Path tempDir;

	@Test
	void fragPipeParameterFileKeepsConditionWithItsAcquisition() throws Exception {
		DDAMagFragPipeParameters params = new DDAMagFragPipeParameters(
				"DDA_MAG_FRAGPIPE",
				new MetaData(new String[] { "C:/data/a.raw", "C:/data/b.raw" }, new String[] { "SampleA", "SampleB" },
						new int[] { 1, 1 }, new int[] { 1, 1 }),
				tempDir.toString(), "", false, "", "HCD-FTMS", 1);
		params.setMetadataTableData(List.of(
				new String[] { "SampleA", "", "", "", "Treated", "", "", "", "", "" },
				new String[] { "SampleB", "", "", "", "Control", "", "", "", "", "" }));

		DDAMagFragPipeParameterIO io = new DDAMagFragPipeParameterIO();
		assertTrue(io.exportParametersToFile(params, tempDir.resolve("params.json").toFile(), "test"));

		DDAMagFragPipeParameters loaded = io.loadParameters(tempDir.resolve("params.json").toFile());
		assertEquals("Treated", loaded.getMetadata().getMetaInfo()[0][0]);
		assertEquals("Control", loaded.getMetadata().getMetaInfo()[1][0]);
		assertEquals("Condition", loaded.getMetadata().getMetaTypeTitles()[0]);

		String report = Files.readString(SampleMetadata.reportFile(tempDir.toFile(), loaded.getMetadata()).toPath());
		assertTrue(report.startsWith("source name\tcomment[data file]\tcomment[fraction identifier]\t"
				+ "comment[technical replicate]\tcharacteristics[condition]\tfactor value[condition]"));
		assertTrue(report.contains("SampleA\ta\t1\t1\tTreated\tTreated"));
		assertTrue(report.contains("SampleB\tb\t1\t1\tControl\tControl"));
	}

	@Test
	void olderPositionalMetadataStillLoadsCondition() {
		JSONObject json = new JSONObject("{\"metadata\":["
				+ "{\"sample_index\":0,\"age\":\"30\",\"gender\":\"F\",\"condition\":\"Treated\"}]}");
		MetaData metadata = new MetaData(new String[] { "a.raw" }, new String[] { "SampleA" }, new int[] { 1 },
				new int[] { 1 });
		SampleMetadata.apply(metadata, SampleMetadata.fromJson(json));
		assertEquals(List.of("Age", "Gender", "Condition"), List.of(metadata.getMetaTypeTitles()));
		assertEquals("Treated", metadata.getMetaInfo()[0][2]);
	}
}
