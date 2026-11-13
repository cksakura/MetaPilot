package uk.ac.quadram.metapilot.resources;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertNull;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Arrays;

import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;

import uk.ac.quadram.metapilot.resources.MetaResourceParameters.ToolType;

class MetaResourceParametersTest {

	@TempDir
	Path tempDir;

	@Test
	void msfraggerIsNotASeparatelyConfigurableDdaResource() {
		assertFalse(Arrays.asList(ToolType.getDDATools()).contains(ToolType.MSFRAGGER));
	}

	@Test
	void resolvesHighestMsfraggerVersionFromFragpipeLauncher() throws Exception {
		Path root = tempDir.resolve("FragPipe-24.0");
		Path launcher = Files.createDirectories(root.resolve("bin")).resolve("fragpipe.bat");
		Files.createFile(launcher);
		createMsfraggerJar(root, "4.4.1");
		Path newestJar = createMsfraggerJar(root, "4.10.0");

		MetaResourceParameters params = new MetaResourceParameters(tempDir.resolve("resources").toString());
		params.setToolPath(ToolType.FRAGPIPE, launcher.toString());

		assertEquals(newestJar.toFile().getCanonicalFile(), params.resolveMsfraggerJar().getCanonicalFile());
		assertTrue(params.isToolAvailable(ToolType.MSFRAGGER));
	}

	@Test
	void missingBundledJarDoesNotSatisfyMsfraggerDependency() throws Exception {
		Path root = tempDir.resolve("FragPipe-24.0");
		Path launcher = Files.createDirectories(root.resolve("bin")).resolve("fragpipe.bat");
		Files.createFile(launcher);

		MetaResourceParameters params = new MetaResourceParameters(tempDir.resolve("resources").toString());
		params.setToolPath(ToolType.FRAGPIPE, launcher.toString());

		assertNull(params.resolveMsfraggerJar());
		assertFalse(params.isToolAvailable(ToolType.MSFRAGGER));
	}

	@Test
	void resolvesVersionedLinuxDiannBinaryForMsbooster() throws Exception {
		Path root = tempDir.resolve("FragPipe-24.0");
		Path launcher = Files.createDirectories(root.resolve("bin")).resolve("fragpipe");
		Files.createFile(launcher);
		Path linux = Files.createDirectories(root.resolve("tools").resolve("diann").resolve("1.8.2_beta_8").resolve("linux"));
		Files.createFile(linux.resolve("libtorch.so"));
		Files.createFile(linux.resolve("unimod.obo"));
		Path binary = Files.createFile(linux.resolve("diann-1.8.1.8"));

		MetaResourceParameters params = new MetaResourceParameters(tempDir.resolve("resources").toString());
		params.setToolPath(ToolType.FRAGPIPE, launcher.toString());

		assertEquals(binary.toFile().getCanonicalFile(), params.resolveMsboosterDiann().getCanonicalFile());
	}

	@Test
	void prefersPlainDiannOverVersionedName() throws Exception {
		Path root = tempDir.resolve("FragPipe-24.0");
		Path launcher = Files.createDirectories(root.resolve("bin")).resolve("fragpipe");
		Files.createFile(launcher);
		Path linux = Files.createDirectories(root.resolve("tools").resolve("diann").resolve("1.8.2_beta_8").resolve("linux"));
		Files.createFile(linux.resolve("diann-1.8.1.8"));
		Path plain = Files.createFile(linux.resolve("diann"));

		MetaResourceParameters params = new MetaResourceParameters(tempDir.resolve("resources").toString());
		params.setToolPath(ToolType.FRAGPIPE, launcher.toString());

		assertEquals(plain.toFile().getCanonicalFile(), params.resolveMsboosterDiann().getCanonicalFile());
	}

	@Test
	void diannBinaryNameMatcherAcceptsFragpipeLinuxLayout() {
		assertTrue(MetaResourceParameters.isMsboosterDiannBinaryName("diann-1.8.1.8"));
		assertTrue(MetaResourceParameters.isMsboosterDiannBinaryName("diann"));
		assertTrue(MetaResourceParameters.isMsboosterDiannBinaryName("diann.exe"));
		assertTrue(MetaResourceParameters.isMsboosterDiannBinaryName("diann-linux"));
		assertFalse(MetaResourceParameters.isMsboosterDiannBinaryName("libtorch.so"));
		assertFalse(MetaResourceParameters.isMsboosterDiannBinaryName("unimod.obo"));
	}

	private static Path createMsfraggerJar(Path fragpipeRoot, String version) throws Exception {
		Path directory = Files.createDirectories(fragpipeRoot.resolve("tools").resolve("MSFragger-" + version));
		return Files.createFile(directory.resolve("MSFragger-" + version + ".jar"));
	}
}
