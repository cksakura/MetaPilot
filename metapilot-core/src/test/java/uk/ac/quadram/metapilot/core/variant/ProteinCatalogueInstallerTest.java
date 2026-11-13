package uk.ac.quadram.metapilot.core.variant;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertNotNull;
import static org.junit.jupiter.api.Assertions.assertNull;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;

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

	@Test
	@DisplayName("The protein_catalogue listing selects the 100% archive and ignores the other cut-offs")
	void listingDetectsTheFullArchiveName() throws IOException {
		String v2 = "http://ftp.ebi.ac.uk/pub/databases/metagenomics/mgnify_genomes/human-gut/v2.0.2";
		List<String> v2names = ProteinCatalogueInstaller.hrefNames(
				"<a href=\"?C=N;O=D\">Name</a> <a href=\"/parent/\">Parent</a>"
						+ "<a href=\"uhgp-50.tar.gz\">uhgp-50.tar.gz</a>"
						+ "<a href=\"uhgp-90.tar.gz\">uhgp-90.tar.gz</a>"
						+ "<a href=\"uhgp-95.tar.gz\">uhgp-95.tar.gz</a>"
						+ "<a href=\"uhgp-100.tar.gz\">uhgp-100.tar.gz</a>");
		assertEquals(v2 + "/protein_catalogue/uhgp-100.tar.gz",
				ProteinCatalogueInstaller.archiveUrl(v2 + "/", v2names));

		String v3 = "http://ftp.ebi.ac.uk/pub/databases/metagenomics/mgnify_genomes/human-gut/v3.0";
		List<String> v3names = ProteinCatalogueInstaller.hrefNames(
				"<a href=\"protein_catalogue-50.tar.gz\">50</a>"
						+ "<a href=\"protein_catalogue-90.tar.gz\">90</a>"
						+ "<a href=\"protein_catalogue-95.tar.gz\">95</a>"
						+ "<a href=\"protein_catalogue-100.tar.gz\">100</a>"
						+ "<a href=\"proteins.faa.gz\">proteins</a>");
		assertEquals(v3 + "/protein_catalogue/protein_catalogue-100.tar.gz",
				ProteinCatalogueInstaller.archiveUrl(v3, v3names));
		assertEquals("renamed-100.tar.gz",
				PangenomeResource.selectFullCatalogueArchive(List.of("renamed-50.tar.gz", "renamed-100.tar.gz",
						"proteins.faa.gz")));
		assertNull(PangenomeResource.selectFullCatalogueArchive(List.of("uhgp-50.tar.gz", "proteins.faa.gz")));
	}

	@Test
	@DisplayName("v2.0.2 uses uhgp-100.tar.gz; v3.0 falls back to protein_catalogue-100.tar.gz")
	void archiveUrlPrefersUhgpThenV3Name() throws IOException {
		String v2 = "http://ftp.ebi.ac.uk/pub/databases/metagenomics/mgnify_genomes/human-gut/v2.0.2";
		assertEquals(v2 + "/protein_catalogue/uhgp-100.tar.gz",
				ProteinCatalogueInstaller.archiveUrl(v2, u -> u.endsWith("uhgp-100.tar.gz")));
		assertEquals(v2 + "/protein_catalogue/uhgp-100.tar.gz",
				ProteinCatalogueInstaller.archiveUrl(v2 + "/", u -> u.endsWith("uhgp-100.tar.gz")));
		String v3 = "http://ftp.ebi.ac.uk/pub/databases/metagenomics/mgnify_genomes/human-gut/v3.0";
		assertEquals(v3 + "/protein_catalogue/protein_catalogue-100.tar.gz",
				ProteinCatalogueInstaller.archiveUrl(v3, u -> u.endsWith("protein_catalogue-100.tar.gz")));
		IOException missing = assertThrows(IOException.class,
				() -> ProteinCatalogueInstaller.archiveUrl(v2, u -> false));
		assertTrue(missing.getMessage().contains("uhgp-100.tar.gz"));
		assertTrue(missing.getMessage().contains("protein_catalogue-100.tar.gz"));
	}

	@Test
	@DisplayName("A renamed 100% archive already on disk is prepared without opening a URL")
	void installFromV3ArchiveName() throws IOException {
		File versionDir = new File(new File(tmp, "human-gut"), "v3.0");
		File dir = ProteinCatalogueInstaller.directory(versionDir);
		dir.mkdirs();
		PangenomeResourcePreparerTest helper = new PangenomeResourcePreparerTest();
		helper.metadata(versionDir);
		helper.archive(dir);
		File v2name = new File(dir, PangenomeResource.ARCHIVE);
		File v3name = new File(dir, "catalogue-release-100.tar.gz");
		assertTrue(v2name.renameTo(v3name));
		assertNotNull(PangenomeResource.locate(null, versionDir));
		assertTrue(PangenomeResource.locate(null, versionDir).needsPrepare());

		assertEquals(dir, ProteinCatalogueInstaller.install(UNREACHABLE, versionDir, true, null));
		assertTrue(ProteinCatalogueInstaller.isInstalled(versionDir));
		assertFalse(v3name.exists());
		assertFalse(v2name.exists());
	}
}
