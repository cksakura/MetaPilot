package uk.ac.quadram.metapilot.task.dda.mag.par;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.io.File;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.List;

import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;

import uk.ac.quadram.metapilot.core.fragpipe.FragpipeModifications;
import uk.ac.quadram.metapilot.core.msfragger.MSFraggerHeadlessRunner;

class DDAMagFragPipeParameterIOTest {

	@TempDir
	Path tempDir;

	private DDAMagFragPipeParameters load(String json) throws Exception {
		File file = tempDir.resolve("params.json").toFile();
		Files.writeString(file.toPath(), json);
		return new DDAMagFragPipeParameterIO().loadParameters(file);
	}

	@Test
	@DisplayName("New parameters search the default modifications with a fully specific digest")
	void defaults() {
		DDAMagFragPipeParameters params = new DDAMagFragPipeParameters();
		assertTrue(params.isApplyModifications());
		assertEquals(FragpipeModifications.DEFAULT_FIXED, params.getFixedModifications());
		assertEquals(FragpipeModifications.DEFAULT_VARIABLE, params.getVariableModifications());
		assertEquals(MSFraggerHeadlessRunner.DIGEST_SPECIFIC, params.getDigestModeString());
	}

	@Test
	@DisplayName("A file written with the marker applies exactly its listed modifications")
	void markedFileAppliesItsModifications() throws Exception {
		DDAMagFragPipeParameters params = load("{\"" + DDAMagFragPipeParameterIO.MODIFICATIONS_APPLIED_KEY
				+ "\": true, \"fixed_modifications\": [\"C:57.02146\"],"
				+ " \"variable_modifications\": [\"Phospho:STY:79.96633\"]}");
		assertTrue(params.isApplyModifications());
		assertEquals(List.of("C:57.02146"), params.getFixedModifications());
		assertEquals(List.of("Phospho:STY:79.96633"), params.getVariableModifications());
	}

	@Test
	@DisplayName("A file from before the marker keeps searching the workflow template's modifications")
	void legacyFileKeepsTemplateModifications() throws Exception {
		DDAMagFragPipeParameters params = load(
				"{\"fixed_modifications\": [\"C:57.02146\"], \"variable_modifications\": []}");
		assertFalse(params.isApplyModifications());
	}
}
