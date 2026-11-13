package uk.ac.quadram.metapilot.core.opensearch;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;

import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;

import uk.ac.quadram.metapilot.core.msfragger.MSFraggerHeadlessRunner;

class OpenSearchSpectraTest {

	@TempDir
	File tmp;

	@Test
	void mzmlInputsAreReadDirectly() throws Exception {
		File data = new File(tmp, "data");
		File a = new File(data, "A.mzML");
		File b = new File(data, "B.mzml");
		File root = new File(tmp, "opensearch_MAG_search");

		assertFalse(OpenSearchSpectra.needsConversion(List.of(a, b)));
		assertEquals(a, OpenSearchSpectra.readable(a, root));
		assertEquals(data, OpenSearchSpectra.spectraDir(List.of(a, b), root));
	}

	@Test
	void vendorInputsAreReadThroughTheirConvertedMzml() throws Exception {
		File d = new File(tmp, "run1.d");
		File raw = new File(tmp, "run2.raw");
		File root = new File(tmp, "opensearch_MAG_search");
		File converted = new File(root, OpenSearchSpectra.CONVERTED_DIR);

		assertTrue(OpenSearchSpectra.needsConversion(List.of(d, raw)));
		assertEquals(new File(converted, "run1.mzML"), OpenSearchSpectra.readable(d, root));
		assertEquals(new File(converted, "run2.mzML"), OpenSearchSpectra.readable(raw, root));
		assertEquals(converted, OpenSearchSpectra.spectraDir(List.of(d, raw), root));
		assertFalse(OpenSearchSpectra.isReadable(d, root), "not converted yet");
	}

	@Test
	void mixedInputsAreRefused() {
		assertThrows(IOException.class, () -> OpenSearchSpectra
				.needsConversion(List.of(new File(tmp, "A.mzML"), new File(tmp, "B.d"))));
	}

	@Test
	void collectMovesMsfraggersUncalibratedMzmlUnderTheRunName() throws Exception {
		File d = new File(tmp, "N2022_lib19.d");
		assertTrue(d.mkdirs());
		File produced = new File(tmp, "N2022_lib19_uncalibrated.mzML");
		Files.writeString(produced.toPath(), "<indexedmzML/>");
		File root = new File(tmp, "opensearch_MAG_search");

		File target = OpenSearchSpectra.collect(d, root);

		assertEquals(new File(new File(root, OpenSearchSpectra.CONVERTED_DIR), "N2022_lib19.mzML"), target);
		assertEquals("<indexedmzML/>", Files.readString(target.toPath()));
		assertFalse(produced.exists(), "moved, not copied: the data drive keeps no second copy");
		assertTrue(OpenSearchSpectra.isReadable(d, root));
		assertFalse(new File(target.getParentFile(), target.getName() + ".part").exists());
	}

	@Test
	void collectFailsWhenMsfraggerWroteNoMzml() {
		File d = new File(tmp, "missing.d");
		assertThrows(IOException.class, () -> OpenSearchSpectra.collect(d, new File(tmp, "root")));
	}

	@Test
	void writeUncalibratedMzmlOnlyChangesTheParamsWhenAsked() throws Exception {
		File fasta = new File(tmp, "db.fasta");
		Files.writeString(fasta.toPath(), ">P1\nPEPTIDEK\n");
		MSFraggerHeadlessRunner.OpenSearchSettings settings = new MSFraggerHeadlessRunner.OpenSearchSettings();
		settings.targetDecoyFasta = fasta;

		String mzml = MSFraggerHeadlessRunner.previewOpenSearchParamsText(settings);
		settings.writeUncalibratedMzml = true;
		String vendor = MSFraggerHeadlessRunner.previewOpenSearchParamsText(settings);

		assertTrue(mzml.matches("(?s).*\\bwrite_uncalibrated_mzml = 0\\b.*"), "mzML runs keep the template value");
		assertTrue(vendor.matches("(?s).*\\bwrite_uncalibrated_mzml = 1\\b.*"));
	}
}
