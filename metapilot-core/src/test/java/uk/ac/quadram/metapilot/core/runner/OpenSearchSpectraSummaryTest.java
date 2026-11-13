package uk.ac.quadram.metapilot.core.runner;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.io.File;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;

import uk.ac.quadram.metapilot.core.mag.MagResultLayoutUtil;
import uk.ac.quadram.metapilot.core.opensearch.OpenSearchLayout;
import uk.ac.quadram.metapilot.core.workflow.ProgressCallback;

/** final_summary.tsv's MS/MS count for runs started from mzML (no conversion folder). */
class OpenSearchSpectraSummaryTest {

	@TempDir
	File tmp;

	private final List<String> messages = new ArrayList<>();
	private final ProgressCallback cb = new ProgressCallback() {
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

	/** An mzML in msconvert's line layout: one MS1 and {@code ms2} MS2 spectra. */
	private File mzml(File dir, String name, int ms2) throws Exception {
		StringBuilder xml = new StringBuilder("<mzML>\n<run>\n<spectrumList>\n");
		int scan = 1;
		xml.append("<spectrum index=\"0\" id=\"controllerType=0 controllerNumber=1 scan=").append(scan++)
				.append("\" defaultArrayLength=\"0\">\n")
				.append("<cvParam cvRef=\"MS\" accession=\"MS:1000579\" name=\"MS1 spectrum\" value=\"\"/>\n")
				.append("<cvParam cvRef=\"MS\" accession=\"MS:1000511\" name=\"ms level\" value=\"1\"/>\n")
				.append("<cvParam cvRef=\"MS\" accession=\"MS:1000016\" name=\"scan start time\" value=\"0.01\"/>\n")
				.append("</spectrum>\n");
		for (int i = 0; i < ms2; i++) {
			xml.append("<spectrum index=\"").append(i + 1).append("\" id=\"controllerType=0 controllerNumber=1 scan=")
					.append(scan++).append("\" defaultArrayLength=\"0\">\n")
					.append("<cvParam cvRef=\"MS\" accession=\"MS:1000580\" name=\"MSn spectrum\" value=\"\"/>\n")
					.append("<cvParam cvRef=\"MS\" accession=\"MS:1000511\" name=\"ms level\" value=\"2\"/>\n")
					.append("<cvParam cvRef=\"MS\" accession=\"MS:1000016\" name=\"scan start time\" value=\"")
					.append(0.02 + i * 0.01).append("\"/>\n</spectrum>\n");
		}
		xml.append("</spectrumList>\n</run>\n</mzML>\n");
		File f = new File(dir, name);
		Files.write(f.toPath(), xml.toString().getBytes(StandardCharsets.UTF_8));
		return f;
	}

	private static List<String[]> experiments(File... mzml) {
		List<String[]> out = new ArrayList<>();
		for (File f : mzml) {
			out.add(new String[] { f.getAbsolutePath(), f.getName().replace(".mzML", "") });
		}
		return out;
	}

	@Test
	void writesMs2CountUnderTheResultRootFromTheMzml() throws Exception {
		File inputs = new File(tmp, "mzML");
		inputs.mkdirs();
		File a = mzml(inputs, "S1.mzML", 3);
		File b = mzml(inputs, "S2.mzML", 5);
		File resultRoot = new File(tmp, "result");

		File dir = MagOpenSearchHeadlessRunner.spectraSummaryDir(experiments(a, b), resultRoot, cb);

		assertEquals(MagResultLayoutUtil.spectraDir(resultRoot), dir);
		assertEquals(Arrays.asList("FileName\tMS2 count", "S1\t3", "S2\t5"),
				Files.readAllLines(new File(dir, "ms2Count.txt").toPath()));
		List<String> rt = Files.readAllLines(new File(dir, "rt.txt").toPath());
		assertEquals(1 + 3 + 5, rt.size(), "one row per MS2 scan, MS1 left out");
		assertEquals("S1\t2\t0.02", rt.get(1));
		assertFalse(new File(inputs, "ms2Count.txt").exists(), "nothing is written beside the inputs");
	}

	@Test
	void finalSummaryUsesTheMs2Count() throws Exception {
		File inputs = new File(tmp, "mzML");
		inputs.mkdirs();
		File a = mzml(inputs, "S1.mzML", 4);
		File dir = MagOpenSearchHeadlessRunner.spectraSummaryDir(experiments(a), new File(tmp, "result"), cb);
		File combined = new File(tmp, "combined");
		File out = new File(tmp, "out");
		java.util.HashMap<String, Integer> fileIdMap = new java.util.HashMap<>();
		fileIdMap.put("S1", 0);
		@SuppressWarnings("unchecked")
		java.util.HashSet<String>[] sequences = new java.util.HashSet[] {
				new java.util.HashSet<>(Arrays.asList("PEPTIDEK")) };
		DdaMagReportExportService.writeIdentificationSummary(combined, out, dir, fileIdMap, new int[] { 2 }, sequences,
				null);
		List<String> summary = Files.readAllLines(new File(out, "final_summary.tsv").toPath());
		assertTrue(summary.get(1).startsWith("S1\t"), summary.get(1));
		String[] row = summary.get(1).split("\t");
		assertEquals("4", row[2], "MS/MS is the MS2 scan count");
		assertEquals("2", row[3], "MS/MS Identified");
		assertEquals("50", row[4], "identification rate");
		assertEquals("2", row[6], "one identification per spectrum: PSMs equal identified spectra");
	}

	@Test
	void identificationRateCountsSpectraAndPsmsAreReportedApart() throws Exception {
		File inputs = new File(tmp, "mzML");
		inputs.mkdirs();
		File a = mzml(inputs, "S1.mzML", 10);
		File dir = MagOpenSearchHeadlessRunner.spectraSummaryDir(experiments(a), new File(tmp, "result"), cb);
		File combined = new File(tmp, "combined");
		File out = new File(tmp, "out");
		java.util.HashMap<String, Integer> fileIdMap = new java.util.HashMap<>();
		fileIdMap.put("S1", 0);
		@SuppressWarnings("unchecked")
		java.util.HashSet<String>[] sequences = new java.util.HashSet[] {
				new java.util.HashSet<>(Arrays.asList("PEPTIDEK", "SAMPLEK")) };
		// 6 PSMs on 5 spectra: one scan carries a second identification (a reassigned precursor).
		DdaMagReportExportService.writeIdentificationSummary(combined, out, dir, fileIdMap, new int[] { 6 },
				new int[] { 5 }, sequences, null);

		assertEquals(Arrays.asList("Raw file\tUnique peptide count\tPSM count\tSpectrum count", "S1\t2\t6\t5"),
				Files.readAllLines(new File(combined, "id.tsv").toPath()));
		List<String> summary = Files.readAllLines(new File(out, "final_summary.tsv").toPath());
		assertEquals("Raw file\tExperiment\tMS/MS\tMS/MS Identified\tMS/MS Identified [%]\t"
				+ "Peptide Sequences Identified\tPSMs", summary.get(0));
		assertEquals(Arrays.asList("S1", "S1", "10", "5", "50", "2", "6"), Arrays.asList(summary.get(1).split("\t")));
		assertEquals(Arrays.asList("Total", "", "10", "5", "50", "2", "6"), Arrays.asList(summary.get(2).split("\t")));

		// Without ms2Count.txt the MS/MS count falls back to the identified spectra, not the PSMs.
		assertTrue(new File(dir, "ms2Count.txt").delete());
		File out2 = new File(tmp, "out2");
		DdaMagReportExportService.writeIdentificationSummary(new File(tmp, "combined2"), out2, null, fileIdMap,
				new int[] { 6 }, new int[] { 5 }, sequences, null);
		assertEquals(Arrays.asList("S1", "S1", "5", "5", "100", "2", "6"),
				Arrays.asList(Files.readAllLines(new File(out2, "final_summary.tsv").toPath()).get(1).split("\t")));
	}

	@Test
	void identifiedSpectraCountEachScanOnce() {
		java.util.HashMap<String, Integer> fileIdMap = new java.util.HashMap<>();
		fileIdMap.put("S1", 0);
		fileIdMap.put("S2", 1);
		List<uk.ac.quadram.metapilot.core.opensearch.OpenSearchPsm> psms = new ArrayList<>();
		for (Object[] f : new Object[][] { { "S1", 7 }, { "S1", 7 }, { "S1", 8 }, { "S2", 7 }, { "S3", 1 } }) {
			uk.ac.quadram.metapilot.core.opensearch.OpenSearchPsm psm = new uk.ac.quadram.metapilot.core.opensearch.OpenSearchPsm();
			psm.setFileName((String) f[0]);
			psm.setScan((Integer) f[1]);
			psms.add(psm);
		}
		int[] spectra = OpenSearchPostAnalysisTsvExporter.identifiedSpectraPerFile(psms, fileIdMap,
				Arrays.asList("S1", "S2"));
		assertEquals(2, spectra[0], "scan 7 of S1 counts once despite two identifications");
		assertEquals(1, spectra[1]);
	}

	@Test
	void prefersAnExistingCountBesideTheInputs() throws Exception {
		File inputs = new File(tmp, "mzML");
		inputs.mkdirs();
		File a = mzml(inputs, "S1.mzML", 3);
		Files.write(new File(inputs, "ms2Count.txt").toPath(), Arrays.asList("FileName\tMS2 count", "S1\t99"));
		File resultRoot = new File(tmp, "result");

		assertEquals(inputs, MagOpenSearchHeadlessRunner.spectraSummaryDir(experiments(a), resultRoot, cb));
		assertFalse(new File(resultRoot, MagOpenSearchHeadlessRunner.SPECTRA_SUMMARY_DIR).exists());
	}

	@Test
	void missingSpectraFallBackWithoutLeavingAPartialFile() throws Exception {
		File inputs = new File(tmp, "mzML");
		inputs.mkdirs();
		File a = mzml(inputs, "S1.mzML", 3);
		File missing = new File(inputs, "S2.mzML");
		File resultRoot = new File(tmp, "result");

		File dir = MagOpenSearchHeadlessRunner.spectraSummaryDir(experiments(a, missing), resultRoot, cb);

		assertEquals(inputs, dir);
		assertFalse(new File(MagResultLayoutUtil.spectraDir(resultRoot), "ms2Count.txt").exists(),
				"a partial count would be taken as complete on the next run");
		assertTrue(messages.stream().anyMatch(m -> m.contains("missing or empty")), messages.toString());
	}

	@Test
	void countsSurviveTheNestedLayoutPassThatUsedToHideThem() throws Exception {
		File inputs = new File(tmp, "mzML");
		inputs.mkdirs();
		File a = mzml(inputs, "S1.mzML", 4);
		File resultRoot = new File(tmp, "opensearch_MAG_search");
		File dir = MagOpenSearchHeadlessRunner.spectraSummaryDir(experiments(a), resultRoot, cb);
		OpenSearchLayout.sampleDbParams(resultRoot, "S1");

		assertEquals(MagResultLayoutUtil.spectraDir(resultRoot), dir);
		assertTrue(new File(dir, "ms2Count.txt").isFile());
		assertFalse(new File(new File(new File(resultRoot, MagResultLayoutUtil.SEARCH_RESULTS_DIR_NAME),
				MagResultLayoutUtil.SPECTRA_SUMMARY_DIR_NAME), "ms2Count.txt").exists());

		File modeRoot = OpenSearchLayout.modeRoot(resultRoot, null);
		writeSummary(modeRoot, dir, 2);
		assertEquals("50", summaryRate(modeRoot));
	}

	@Test
	void restoresACopyAlreadyNestedUnderSearchResults() throws Exception {
		File inputs = new File(tmp, "mzML");
		inputs.mkdirs();
		File a = mzml(inputs, "S1.mzML", 4);
		File resultRoot = new File(tmp, "opensearch_MAG_search");
		File nested = new File(new File(resultRoot, MagResultLayoutUtil.SEARCH_RESULTS_DIR_NAME),
				MagResultLayoutUtil.SPECTRA_SUMMARY_DIR_NAME);
		assertTrue(nested.mkdirs());
		Files.write(new File(nested, "ms2Count.txt").toPath(), Arrays.asList("FileName\tMS2 count", "S1\t4"));

		File dir = MagOpenSearchHeadlessRunner.spectraSummaryDir(experiments(a), resultRoot, cb);

		assertEquals(MagResultLayoutUtil.spectraDir(resultRoot), dir);
		assertTrue(new File(dir, "ms2Count.txt").isFile());
		OpenSearchLayout.sampleDbParams(resultRoot, "S1");
		assertTrue(new File(dir, "ms2Count.txt").isFile());
	}

	@Test
	void identificationRateFindsCountsInTheResultSpectraFolder() throws Exception {
		File engine = new File(tmp, "opensearch_MAG_search");
		File spectra = new File(tmp, "spectra");
		assertTrue(spectra.mkdirs());
		Files.write(new File(spectra, "ms2Count.txt").toPath(), Arrays.asList("FileName\tMS2 count", "S1\t4"));
		Files.write(new File(spectra, "rt.txt").toPath(), Arrays.asList("FileName\tScannum\tRt", "S1\t2\t0.02"));

		assertEquals(new File(spectra, "ms2Count.txt"),
				DdaMagReportExportService.resolveMs2CountFile(null, engine));
		assertEquals(new File(spectra, "ms2Count.txt"),
				DdaMagReportExportService.resolveMs2CountFile(new File(tmp, "mzML"), engine));
		writeSummary(engine, null, 2);
		assertEquals("50", summaryRate(engine));
	}

	@Test
	void identificationRateFindsALeftoverSpectraSummary() throws Exception {
		File engine = new File(tmp, "fragpipe_MAG_search");
		File legacy = new File(engine, MagResultLayoutUtil.SPECTRA_SUMMARY_DIR_NAME);
		assertTrue(legacy.mkdirs());
		Files.write(new File(legacy, "ms2Count.txt").toPath(), Arrays.asList("FileName\tMS2 count", "S1\t8"));

		assertEquals(new File(legacy, "ms2Count.txt"),
				DdaMagReportExportService.resolveMs2CountFile(null, engine));
		writeSummary(engine, null, 2);
		assertEquals("25", summaryRate(engine));
	}

	@Test
	void writerFindsCountsAfterTheyWereNested() throws Exception {
		File engine = new File(tmp, "opensearch_MAG_search");
		File nested = new File(new File(engine, MagResultLayoutUtil.SEARCH_RESULTS_DIR_NAME),
				MagResultLayoutUtil.SPECTRA_SUMMARY_DIR_NAME);
		assertTrue(nested.mkdirs());
		Files.write(new File(nested, "ms2Count.txt").toPath(), Arrays.asList("FileName\tMS2 count", "S1\t4"));
		File modeRoot = OpenSearchLayout.modeRoot(engine, null);
		writeSummary(modeRoot, new File(engine, MagResultLayoutUtil.SPECTRA_SUMMARY_DIR_NAME), 2);
		assertEquals("50", summaryRate(modeRoot));
	}

	@Test
	void fragPipeHpcWritesCountsUnderTheEngineRootWhenThereIsNoConversionFolder() throws Exception {
		File inputs = new File(tmp, "mzML");
		inputs.mkdirs();
		File a = mzml(inputs, "S1.mzML", 4);
		File resultRoot = new File(tmp, "fragpipe_MAG_search");
		File dir = MagFragPipeHeadlessRunner.prepareSpectraSummaryDir(resultRoot, null,
				new String[] { a.getAbsolutePath() }, cb);
		assertEquals(MagResultLayoutUtil.spectraDir(resultRoot), dir);
		assertTrue(new File(dir, "ms2Count.txt").isFile());
		MagResultLayoutUtil.resolveSearchResultsRoot(resultRoot);
		assertTrue(new File(dir, "ms2Count.txt").isFile());
	}

	private void writeSummary(File modeRoot, File spectraDir, int identified) throws Exception {
		File combined = new File(modeRoot, "combined");
		java.util.HashMap<String, Integer> fileIdMap = new java.util.HashMap<>();
		fileIdMap.put("S1", 0);
		@SuppressWarnings("unchecked")
		java.util.HashSet<String>[] sequences = new java.util.HashSet[] {
				new java.util.HashSet<>(Arrays.asList("PEPTIDEK")) };
		DdaMagReportExportService.writeIdentificationSummary(combined, modeRoot, spectraDir, fileIdMap,
				new int[] { identified }, sequences, null);
	}

	@Test
	void fragpipeRenamedRunsFindTheirMs2Count() {
		java.util.Map<String, Integer> counts = java.util.Map.of("Qiao-FB1-DDA-1", 85389, "Qiao-FB1-DDA-2", 85471,
				"YS-194-1-1", 142931);
		assertEquals(85389, DdaMagReportExportService.ms2CountOf(counts, "Qiao-FB1-DDA-1"), "spectrum-file name");
		assertEquals(85471, DdaMagReportExportService.ms2CountOf(counts, "Qiao_FB1_DDA_2_1"), "FragPipe name");
		assertEquals(142931, DdaMagReportExportService.ms2CountOf(counts, "YS_194_1_1_1"));
		assertEquals(null, DdaMagReportExportService.ms2CountOf(counts, "Qiao_FB9_DDA_1_1"));
	}

	private String summaryRate(File modeRoot) throws Exception {
		String[] row = Files.readAllLines(new File(modeRoot, "final_summary.tsv").toPath()).get(1).split("\t");
		return row[4];
	}
}
