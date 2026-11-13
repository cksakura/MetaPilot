package uk.ac.quadram.metapilot.core.runner;

import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.io.File;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;

import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;

import uk.ac.quadram.metapilot.core.opensearch.OpenSearchLayout;
import uk.ac.quadram.metapilot.core.workflow.ProgressCallback;

class MagOpenSearchIndexCleanupTest {

	@TempDir
	Path tmp;

	@Test
	void removesPepindexBesideTheRefinedDatabaseAndLeavesTheFasta() throws Exception {
		File resultRoot = tmp.resolve("result").toFile();
		Files.createDirectories(resultRoot.toPath());
		Path marker = OpenSearchLayout.markerDir(resultRoot, "S1").toPath();
		Files.createDirectories(marker);
		Path fasta = marker.resolve("S1.fasta");
		Path index = marker.resolve("S1.fasta.1.pepindex");
		Path pepxml = marker.resolve("S1.pepXML");
		Files.writeString(fasta, ">p\nPEPTIDEK\n");
		Files.writeString(index, "index");
		Files.writeString(pepxml, "<msms_pipeline_analysis/>");
		Path fallbackIndex = resultRoot.toPath().resolve("marker_target_decoy.fasta.1.pepindex");
		Files.writeString(fallbackIndex, "fallback");
		Path magDbIndex = tmp.resolve("magdb").resolve("Marker4_td.fasta.1.pepindex");
		Files.createDirectories(magDbIndex.getParent());
		Files.writeString(magDbIndex, "shared");

		List<String> messages = new ArrayList<>();
		MagOpenSearchHeadlessRunner.cleanupDatabaseIndexes(resultRoot,
				List.<String[]>of(new String[] { "S1.mzML" }), recording(messages));

		assertFalse(Files.exists(index));
		assertFalse(Files.exists(fallbackIndex));
		assertTrue(Files.exists(fasta));
		assertTrue(Files.exists(pepxml));
		assertTrue(Files.exists(magDbIndex), "shared MAG-database index must not be deleted");
		assertTrue(messages.stream().anyMatch(m -> m.contains("Deleted MSFragger database index")));
	}

	@Test
	void alsoClearsIndexesCopiedIntoSampleDbFolders() throws Exception {
		File resultRoot = tmp.toFile();
		Path sampleDb = OpenSearchLayout.sampleDbDir(resultRoot, "S1").toPath();
		Files.createDirectories(sampleDb);
		Path index = sampleDb.resolve("S1.fasta.2.pepindex");
		Files.writeString(index, "index");

		MagOpenSearchHeadlessRunner.cleanupDatabaseIndexes(resultRoot, List.<String[]>of(new String[] { "S1.mzML" }),
				new ProgressCallback.NoOp());

		assertFalse(Files.exists(index));
	}

	private static ProgressCallback recording(List<String> messages) {
		return new ProgressCallback() {
			@Override
			public void onProgress(double workDone, double max) {
			}

			@Override
			public void onMessage(String message) {
				messages.add(message);
			}

			@Override
			public void onComplete(Object result) {
			}

			@Override
			public void onError(Exception error) {
			}
		};
	}
}
