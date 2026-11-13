package uk.ac.quadram.metapilot.core.msfragger;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.nio.file.Files;
import java.nio.file.Path;

import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;

class MSFraggerHeadlessRunnerTest {

	@TempDir
	Path tempDir;

	@Test
	void usesBundledWindowsJavaExecutable() throws Exception {
		Path javaExe = Files.createDirectories(tempDir.resolve("runtime").resolve("bin")).resolve("java.exe");
		Files.createFile(javaExe);

		assertEquals(javaExe.toAbsolutePath().normalize().toString(),
				MSFraggerHeadlessRunner.resolveJavaCommand(tempDir.resolve("runtime").toString(), "Windows 11"));
	}

	@Test
	void usesEclipseLaunchJvm() {
		String command = MSFraggerHeadlessRunner.resolveJavaCommand();

		assertTrue(Path.of(command).isAbsolute(), "Expected the Java executable from the current JVM");
		assertTrue(Files.isRegularFile(Path.of(command)), "Current JVM Java executable should exist");
	}

	@Test
	void fallsBackToPathWhenJavaHomeHasNoLauncher() {
		assertEquals("java", MSFraggerHeadlessRunner.resolveJavaCommand(tempDir.resolve("missing").toString(),
				"Windows 11"));
	}

	@Test
	void digestionModeSetsEnzymeTermini() {
		assertEquals(2, MSFraggerHeadlessRunner.enzymeTermini(MSFraggerHeadlessRunner.DIGEST_SPECIFIC, "Trypsin"));
		assertEquals(1, MSFraggerHeadlessRunner.enzymeTermini(MSFraggerHeadlessRunner.DIGEST_SEMI_SPECIFIC, "Trypsin"));
		assertEquals(0, MSFraggerHeadlessRunner.enzymeTermini(MSFraggerHeadlessRunner.DIGEST_NONSPECIFIC, "Trypsin"));
		assertEquals(2, MSFraggerHeadlessRunner.enzymeTermini(null, "Trypsin"));
		assertEquals(1, MSFraggerHeadlessRunner.enzymeTermini("Semi free N-term", "Trypsin"));
		assertEquals(0, MSFraggerHeadlessRunner.enzymeTermini(MSFraggerHeadlessRunner.DIGEST_SPECIFIC, "Nonspecific"));
		assertThrows(IllegalArgumentException.class, () -> MSFraggerHeadlessRunner.enzymeTermini("Fully", "Trypsin"));
	}

	@Test
	void aspNCleavesBeforeItsResidue() {
		assertEquals("N", MSFraggerHeadlessRunner.enzymeSense("AspN"));
		assertEquals("N", MSFraggerHeadlessRunner.enzymeSense("Asp-N"));
		assertEquals("C", MSFraggerHeadlessRunner.enzymeSense("Trypsin"));
		assertEquals("C", MSFraggerHeadlessRunner.enzymeSense("GluC"));
	}
}
