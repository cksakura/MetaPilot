package uk.ac.quadram.metapilot.quant.flashLFQ;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.io.File;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.List;

import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;

class FlashLfqTaskTest {

	@TempDir
	Path tempDir;

	@Test
	void detectsWindowsOs() {
		assertTrue(FlashLfqTask.isWindowsOs("Windows 11"));
		assertTrue(FlashLfqTask.isWindowsOs("Windows 10"));
		assertFalse(FlashLfqTask.isWindowsOs("Linux"));
		assertFalse(FlashLfqTask.isWindowsOs("Mac OS X"));
	}

	@Test
	void resolvesDirectoryToCmdExeOnWindows() throws Exception {
		Path dir = Files.createDirectories(tempDir.resolve("flashlfq"));
		Files.createFile(dir.resolve("CMD.exe"));
		Files.createFile(dir.resolve("CMD.dll"));

		File artifact = FlashLfqTask.resolveFlashLfqArtifact(dir.toFile(), true);
		assertEquals("CMD.exe", artifact.getName());
	}

	@Test
	void resolvesDirectoryToCmdDllOnLinux() throws Exception {
		Path dir = Files.createDirectories(tempDir.resolve("flashlfq"));
		Files.createFile(dir.resolve("CMD.exe"));
		Files.createFile(dir.resolve("CMD.dll"));

		File artifact = FlashLfqTask.resolveFlashLfqArtifact(dir.toFile(), false);
		assertEquals("CMD.dll", artifact.getName());
	}

	@Test
	void prefersSiblingDllWhenExeConfiguredOnLinux() throws Exception {
		Path dir = Files.createDirectories(tempDir.resolve("flashlfq"));
		Path exe = Files.createFile(dir.resolve("CMD.exe"));
		Files.createFile(dir.resolve("CMD.dll"));

		File artifact = FlashLfqTask.resolveFlashLfqArtifact(exe.toFile(), false);
		assertEquals("CMD.dll", artifact.getName());
	}

	@Test
	void windowsCommandUsesCmdExe() throws Exception {
		Path exe = Files.createFile(tempDir.resolve("CMD.exe"));
		List<String> cmd = FlashLfqTask.buildFlashLfqCommand(exe.toFile(), true, "/data/psms.tsv", "/data/raw",
				"/data/out", true, true, 8, "dotnet");

		assertEquals("cmd.exe", cmd.get(0));
		assertEquals("/c", cmd.get(1));
		assertEquals("call", cmd.get(2));
		assertTrue(cmd.get(3).endsWith("CMD.exe"));
		assertEquals("--idt", cmd.get(4));
		assertEquals("\"/data/psms.tsv\"", cmd.get(5));
		assertEquals("--mbr", cmd.get(10));
		assertEquals("--nor", cmd.get(11));
		assertEquals("--thr", cmd.get(12));
		assertEquals("8", cmd.get(13));
	}

	@Test
	void linuxCommandUsesDotnetDll() throws Exception {
		Path dll = Files.createFile(tempDir.resolve("CMD.dll"));
		List<String> cmd = FlashLfqTask.buildFlashLfqCommand(dll.toFile(), false, "/data/psms.tsv", "/data/raw",
				"/data/out", true, false, 30, "/opt/dotnet/dotnet");

		assertEquals("/opt/dotnet/dotnet", cmd.get(0));
		assertTrue(cmd.get(1).endsWith("CMD.dll"));
		assertEquals("--idt", cmd.get(2));
		assertEquals("/data/psms.tsv", cmd.get(3));
		assertEquals("--rep", cmd.get(4));
		assertEquals("/data/raw", cmd.get(5));
		assertEquals("--out", cmd.get(6));
		assertEquals("/data/out", cmd.get(7));
		// Switches: present means on, so a false option is left out rather than written "false".
		assertEquals("--mbr", cmd.get(8));
		assertEquals("--thr", cmd.get(9));
		assertEquals("30", cmd.get(10));
		assertFalse(cmd.contains("--nor"));
		assertFalse(cmd.contains("false"));
	}

	@Test
	void windowsDllAlsoUsesDotnet() throws Exception {
		Path dll = Files.createFile(tempDir.resolve("CMD.dll"));
		List<String> cmd = FlashLfqTask.buildFlashLfqCommand(dll.toFile(), true, "C:\\psms.tsv", "C:\\raw", "C:\\out",
				false, true, 4, "dotnet");

		assertEquals("dotnet", cmd.get(0));
		assertTrue(cmd.get(1).endsWith("CMD.dll"));
		assertEquals("C:\\psms.tsv", cmd.get(3));
	}
}
