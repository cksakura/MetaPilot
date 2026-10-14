package uk.ac.quadram.metapilot.core.variant;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertNull;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;

import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;

import uk.ac.quadram.metapilot.core.MGYG.MGYGMagDownloader;

/**
 * The offline paths of the installer: the archive already on disk, so no URL is
 * ever opened ({@code UNREACHABLE} would fail if it were).
 */
class ProteinCatalogueInstallerTest {

	private static final String UNREACHABLE = "http://unreachable.invalid/human-gut/v2.0.2";

	@TempDir
	File tmp;

	private File versionDir() throws IOException {
		File versionDir = new File(new File(tmp, "human-gut"), "v2.0.2");
		File dir = ProteinCatalogueInstaller.directory(versionDir);
		dir.mkdirs();
		PangenomeResourcePreparerTest helper = new PangenomeResourcePreparerTest();
		helper.metadata(versionDir);
		helper.archive(dir);
		return versionDir;
	}

	@Test
	@DisplayName("Human-gut only")
	void supports() {
		assertTrue(ProteinCatalogueInstaller.supports("human-gut"));
		assertTrue(ProteinCatalogueInstaller.supports(" Human-Gut "));
		assertFalse(ProteinCatalogueInstaller.supports("marine"));
		assertFalse(ProteinCatalogueInstaller.supports(null));
	}

	@Test
	@DisplayName("Prepares from an archive on disk, deletes it, and a second call is a no-op")
	void installDeletingArchive() throws IOException {
		File versionDir = versionDir();
		File dir = ProteinCatalogueInstaller.directory(versionDir);
		File archive = new File(dir, PangenomeResource.ARCHIVE);
		assertFalse(ProteinCatalogueInstaller.isInstalled(versionDir));

		assertEquals(dir, ProteinCatalogueInstaller.install(UNREACHABLE, versionDir, true, null));
		assertTrue(ProteinCatalogueInstaller.isInstalled(versionDir));
		assertFalse(archive.exists());
		assertTrue(new File(dir, PangenomeResource.REP_SPECIES).isFile());

		long size = new File(dir, PangenomeResource.FAA).length();
		ProteinCatalogueInstaller.install(UNREACHABLE, versionDir, true, null);
		assertEquals(size, new File(dir, PangenomeResource.FAA).length());
	}

	@Test
	@DisplayName("Keeps the archive when asked")
	void installKeepingArchive() throws IOException {
		File versionDir = versionDir();
		ProteinCatalogueInstaller.install(UNREACHABLE, versionDir, false, null);
		assertTrue(ProteinCatalogueInstaller.isInstalled(versionDir));
		assertTrue(new File(ProteinCatalogueInstaller.directory(versionDir), PangenomeResource.ARCHIVE).isFile());
	}

	@Test
	@DisplayName("A partial download is never mistaken for an archive to prepare")
	void partialDownloadIgnored() throws IOException {
		File versionDir = new File(new File(tmp, "human-gut"), "v2.0.2");
		File dir = ProteinCatalogueInstaller.directory(versionDir);
		dir.mkdirs();
		Files.write(new File(dir, PangenomeResource.ARCHIVE + ".download").toPath(), new byte[100]);
		assertNull(PangenomeResource.locate("", versionDir));
		assertFalse(ProteinCatalogueInstaller.isInstalled(versionDir));
	}

	@Test
	@DisplayName("Preparation-only entry point rejects other catalogues and missing databases")
	void prepareOnlyChecks() {
		assertThrows(IOException.class,
				() -> MGYGMagDownloader.prepareProteinCatalogue("marine", "v1.0", tmp, true, null));
		assertThrows(IOException.class,
				() -> MGYGMagDownloader.prepareProteinCatalogue("human-gut", "v9.9", tmp, true, null));
	}

	@Test
	@DisplayName("Prepared size estimate scales with the archive")
	void preparedEstimate() {
		assertEquals(ProteinCatalogueInstaller.PREPARED_BYTES,
				ProteinCatalogueInstaller.preparedBytes(ProteinCatalogueInstaller.ARCHIVE_BYTES), 1.0);
	}
}
