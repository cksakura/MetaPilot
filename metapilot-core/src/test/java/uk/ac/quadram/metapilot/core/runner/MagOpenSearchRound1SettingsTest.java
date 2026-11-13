package uk.ac.quadram.metapilot.core.runner;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertNotEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.io.File;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.regex.Pattern;

import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;

import uk.ac.quadram.metapilot.core.msfragger.MSFraggerHeadlessRunner;
import uk.ac.quadram.metapilot.task.dda.mag.par.DDAMagOpenSearchParameterIO;
import uk.ac.quadram.metapilot.task.dda.mag.par.DDAMagOpenSearchParameters;

class MagOpenSearchRound1SettingsTest {

	@TempDir
	Path tmp;

	private File fasta() throws Exception {
		return Files.writeString(tmp.resolve("marker_td.fasta"), ">p1\nPEPTIDEK\n").toFile();
	}

	private static boolean hasLine(String text, String key, String value) {
		return Pattern.compile("(?m)^" + Pattern.quote(key) + "\\s*=\\s*" + Pattern.quote(value) + "\\s*(#.*)?$")
				.matcher(text).find();
	}

	@Test
	void closedIsTheDefault() {
		DDAMagOpenSearchParameters params = new DDAMagOpenSearchParameters();
		assertTrue(params.isRound1ClosedSearch());
		assertEquals(DDAMagOpenSearchParameters.ROUND1_SEARCH_CLOSED, params.getRound1SearchMode());
	}

	@Test
	void guiLockedDefaultsRestoreRankDepthAndRound1Mode() {
		DDAMagOpenSearchParameters params = new DDAMagOpenSearchParameters();
		params.setRound1OutputReportTopN(5);
		params.setRound2OutputReportTopN(1);
		params.setClosedStageTopN(1);
		params.setOutputReportTopN(1);
		params.setRound1SearchMode(DDAMagOpenSearchParameters.ROUND1_SEARCH_OPEN);
		params.setClosedStageEnabled(false);
		params.setPrecursorMassLower(-20);
		params.setPrecursorMassUpper(20);
		params.setLocalFdr(0.05);

		params.applyGuiLockedDefaults();

		assertEquals(1, params.getRound1OutputReportTopN());
		assertEquals(3, params.getRound2OutputReportTopN());
		assertEquals(3, params.getClosedStageTopN());
		assertEquals(5, params.getOutputReportTopN());
		assertEquals(DDAMagOpenSearchParameters.ROUND1_SEARCH_CLOSED, params.getRound1SearchMode());
		assertTrue(params.isClosedStageEnabled());
		assertEquals(-150.0, params.getPrecursorMassLower());
		assertEquals(500.0, params.getPrecursorMassUpper());
		assertEquals(0.01, params.getLocalFdr());
		assertEquals(DDAMagOpenSearchParameters.ENGINES_MSFRAGGER, params.getSearchEngines());
	}

	@Test
	void entrapmentIsOffByDefaultAndSurvivesGuiLockedDefaults() {
		DDAMagOpenSearchParameters params = new DDAMagOpenSearchParameters();
		assertFalse(params.isEntrapmentEnabled());
		assertEquals(1L, params.getEntrapmentSeed());
		assertEquals(8, params.getEntrapmentIdDigit());
		params.setEntrapmentEnabled(true);
		params.setEntrapmentSeed(7);
		params.setEntrapmentIdDigit(3);
		params.applyGuiLockedDefaults();
		assertTrue(params.isEntrapmentEnabled());
		assertEquals(7L, params.getEntrapmentSeed());
		assertEquals(3, params.getEntrapmentIdDigit());
	}

	@Test
	void entrapmentFlagRoundTripsThroughParameterJson() throws Exception {
		DDAMagOpenSearchParameters params = new DDAMagOpenSearchParameters();
		params.setResult(tmp.toString());
		params.setEntrapmentEnabled(true);
		params.setEntrapmentSeed(2);
		params.setEntrapmentIdDigit(5);
		File out = tmp.resolve("mag_opensearch_parameters.json").toFile();
		assertTrue(new DDAMagOpenSearchParameterIO().saveParameters(params, out));
		DDAMagOpenSearchParameters loaded = new DDAMagOpenSearchParameterIO().loadParameters(out);
		assertTrue(loaded.isEntrapmentEnabled());
		assertEquals(2L, loaded.getEntrapmentSeed());
		assertEquals(5, loaded.getEntrapmentIdDigit());
	}

	@Test
	void guiLockedDefaultsLeaveGenomeSelectionAndEnzymeAlone() {
		DDAMagOpenSearchParameters params = new DDAMagOpenSearchParameters();
		params.setEnzyme("Trypsin/P");
		params.setMissCleavages(3);
		params.setGenomeSelectorMode(uk.ac.quadram.metapilot.task.mag.MagDbItem.MARKER4);
		params.setGenomeSelectionStrategy(
				uk.ac.quadram.metapilot.task.mag.GenomeSelectionCut.STRATEGY_PSM_COVERAGE);
		params.setGenomeSelectionPsmCoverage(0.99);
		params.setGenomeSelectionStepGuard(0.0);
		params.setGenomeSelectionMaxGenomes(500);
		params.setThreadCount(8);
		params.setMsfraggerMaxHeap("16g");
		params.applyGuiLockedDefaults();
		assertEquals("Trypsin/P", params.getEnzyme());
		assertEquals(3, params.getMissCleavages());
		assertEquals(uk.ac.quadram.metapilot.task.mag.MagDbItem.MARKER4, params.getGenomeSelectorMode());
		assertEquals(uk.ac.quadram.metapilot.task.mag.GenomeSelectionCut.STRATEGY_PSM_COVERAGE,
				params.getGenomeSelectionStrategy());
		assertEquals(0.99, params.getGenomeSelectionPsmCoverage(), 1e-12);
		assertEquals(0.0, params.getGenomeSelectionStepGuard(), 1e-12);
		assertEquals(500, params.getGenomeSelectionMaxGenomes());
		assertEquals(8, params.getThreadCount());
		assertEquals("16g", params.getMsfraggerMaxHeap());

		uk.ac.quadram.metapilot.core.mag.MagGenomeSelectionService.Request request =
				new uk.ac.quadram.metapilot.core.mag.MagGenomeSelectionService.Request();
		MagOpenSearchHeadlessRunner.copyGenomeSelection(request, params);
		assertEquals(uk.ac.quadram.metapilot.task.mag.MagDbItem.MARKER4, request.genomeSelectorMode);
		assertEquals(uk.ac.quadram.metapilot.task.mag.GenomeSelectionCut.STRATEGY_PSM_COVERAGE,
				request.genomeSelectionStrategy);
		assertEquals(0.99, request.genomeSelectionPsmCoverage, 1e-12);
		assertEquals(0.0, request.genomeSelectionStepGuard, 1e-12);
		assertEquals(500, request.genomeSelectionMaxGenomes);
		assertEquals("PSM coverage, target 0.99, step guard 0.0, max genomes 500",
				MagOpenSearchHeadlessRunner.genomeSelectionSummary(params));
	}

	@Test
	void sageEngineModeIsCoercedToMsfragger() {
		DDAMagOpenSearchParameters params = new DDAMagOpenSearchParameters();
		params.setSearchEngines(DDAMagOpenSearchParameters.ENGINES_MSFRAGGER_AND_SAGE);
		assertEquals(DDAMagOpenSearchParameters.ENGINES_MSFRAGGER, params.getSearchEngines());
		assertFalse(params.usesSage());
		assertTrue(params.isMsfraggerOnly());
	}

	@Test
	void closedSearchWritesPpmWindowIsotopeErrorsAndVariableMods() throws Exception {
		DDAMagOpenSearchParameters params = new DDAMagOpenSearchParameters();
		params.setCrystalcEnabled(true);

		MSFraggerHeadlessRunner.OpenSearchSettings settings = MagOpenSearchHeadlessRunner.round1Settings(params,
				fasta());
		String text = MSFraggerHeadlessRunner.previewOpenSearchParamsText(settings);

		assertTrue(hasLine(text, "precursor_mass_units", "1"), text);
		assertTrue(hasLine(text, "precursor_mass_lower", "-20.0"), text);
		assertTrue(hasLine(text, "precursor_mass_upper", "20.0"), text);
		assertTrue(hasLine(text, "isotope_error", "0/1/2"), text);
		assertTrue(hasLine(text, "variable_mod_01", "15.9949 M 3"), text);
		assertTrue(hasLine(text, "variable_mod_02", "42.0106 [^ 1"), text);
		assertFalse(MagOpenSearchHeadlessRunner.round1UsesCrystalC(params));
	}

	@Test
	void openSearchKeepsTheSharedDaWindowAndCrystalC() throws Exception {
		DDAMagOpenSearchParameters params = new DDAMagOpenSearchParameters();
		params.setRound1SearchMode(DDAMagOpenSearchParameters.ROUND1_SEARCH_OPEN);
		params.setCrystalcEnabled(true);

		MSFraggerHeadlessRunner.OpenSearchSettings settings = MagOpenSearchHeadlessRunner.round1Settings(params,
				fasta());
		String text = MSFraggerHeadlessRunner.previewOpenSearchParamsText(settings);

		assertEquals(0, settings.precursorMassUnits);
		assertEquals(params.getPrecursorMassLower(), settings.precursorMassLower);
		assertEquals(params.getPrecursorMassUpper(), settings.precursorMassUpper);
		assertTrue(hasLine(text, "isotope_error", "0"), text);
		assertFalse(text.contains("variable_mod_01"), text);
		assertTrue(MagOpenSearchHeadlessRunner.round1UsesCrystalC(params));
	}

	@Test
	void switchingModeChangesTheParamsText() throws Exception {
		DDAMagOpenSearchParameters closed = new DDAMagOpenSearchParameters();
		DDAMagOpenSearchParameters open = new DDAMagOpenSearchParameters();
		open.setRound1SearchMode(DDAMagOpenSearchParameters.ROUND1_SEARCH_OPEN);
		File fasta = fasta();

		assertNotEquals(
				MSFraggerHeadlessRunner.previewOpenSearchParamsText(MagOpenSearchHeadlessRunner.round1Settings(closed,
						fasta)),
				MSFraggerHeadlessRunner.previewOpenSearchParamsText(MagOpenSearchHeadlessRunner.round1Settings(open,
						fasta)));
	}

	@Test
	void closedStageSearchesTheRefinedDatabaseClosed() throws Exception {
		DDAMagOpenSearchParameters params = new DDAMagOpenSearchParameters();
		assertTrue(params.isClosedStageEnabled());
		// Both searches keep the top 3 candidates by default.
		assertEquals(3, params.getClosedStageTopN());
		assertEquals(3, params.getRound2OutputReportTopN());

		// Every integration mode searches closed at closed_stage_topn.
		assertTrue(params.isJointIdentification());
		assertEquals(3, MagOpenSearchHeadlessRunner.closedStageSettings(params, fasta()).outputReportTopN);
		params.setClosedOpenIntegration(DDAMagOpenSearchParameters.INTEGRATION_CASCADE);
		MSFraggerHeadlessRunner.OpenSearchSettings settings = MagOpenSearchHeadlessRunner.closedStageSettings(params,
				fasta());
		String text = MSFraggerHeadlessRunner.previewOpenSearchParamsText(settings);

		assertEquals(3, settings.outputReportTopN);
		assertTrue(hasLine(text, "output_report_topN", "3"), text);
		assertTrue(hasLine(text, "precursor_mass_units", "1"), text);
		assertTrue(hasLine(text, "precursor_mass_lower", "-20.0"), text);
		assertTrue(hasLine(text, "precursor_mass_upper", "20.0"), text);
		assertTrue(hasLine(text, "isotope_error", "0/1/2"), text);
		assertTrue(hasLine(text, "variable_mod_01", "15.9949 M 3"), text);
	}

	@Test
	void identificationModeStampChangesOnlyWithTheSwitch() throws Exception {
		DDAMagOpenSearchParameters params = new DDAMagOpenSearchParameters();
		File root = tmp.toFile();
		File stamp = MagOpenSearchHeadlessRunner.writeIdentificationModeStamp(params, root, null);
		assertTrue(stamp.isFile());
		assertTrue(stamp.setLastModified(1_000_000L));
		MagOpenSearchHeadlessRunner.writeIdentificationModeStamp(params, root, null);
		assertEquals(1_000_000L, stamp.lastModified());

		params.setClosedStageEnabled(false);
		MagOpenSearchHeadlessRunner.writeIdentificationModeStamp(params, root, null);
		assertNotEquals(1_000_000L, stamp.lastModified());
		assertEquals("open_only", Files.readString(stamp.toPath()).trim());
	}

	@Test
	void identificationModeStampChangesWhenDiannBecomesAvailable() throws Exception {
		DDAMagOpenSearchParameters params = new DDAMagOpenSearchParameters();
		params.setClosedStageEnabled(true);
		params.setClosedOpenIntegration(DDAMagOpenSearchParameters.INTEGRATION_JOINT);
		File root = tmp.resolve("ms2").toFile();
		params.setDiannExecutablePath(new File(root, "missing-diann.exe").getAbsolutePath());
		File stamp = MagOpenSearchHeadlessRunner.writeIdentificationModeStamp(params, root, null);
		String before = Files.readString(stamp.toPath()).trim();
		assertTrue(before.contains("ms2_inputs=no_diann"), before);

		File diann = new File(root, "DiaNN.exe");
		Files.write(diann.toPath(), new byte[] { 0 });
		params.setDiannExecutablePath(diann.getAbsolutePath());
		MagOpenSearchHeadlessRunner.writeIdentificationModeStamp(params, root, null);
		String after = Files.readString(stamp.toPath()).trim();
		assertNotEquals(before, after, "installing DIA-NN must make the joint validation stale");
	}

	@Test
	void unknownModeFallsBackToClosed() {
		DDAMagOpenSearchParameters params = new DDAMagOpenSearchParameters();
		params.setRound1SearchMode("something");
		assertTrue(params.isRound1ClosedSearch());
	}
}
