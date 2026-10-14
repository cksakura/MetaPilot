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
	void closedStageSearchesTheRefinedDatabaseClosedWithOneCandidate() throws Exception {
		DDAMagOpenSearchParameters params = new DDAMagOpenSearchParameters();
		assertTrue(params.isClosedStageEnabled());
		params.setRound2OutputReportTopN(5);

		MSFraggerHeadlessRunner.OpenSearchSettings settings = MagOpenSearchHeadlessRunner.closedStageSettings(params,
				fasta());
		String text = MSFraggerHeadlessRunner.previewOpenSearchParamsText(settings);

		assertEquals(1, settings.outputReportTopN);
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
	void unknownModeFallsBackToClosed() {
		DDAMagOpenSearchParameters params = new DDAMagOpenSearchParameters();
		params.setRound1SearchMode("something");
		assertTrue(params.isRound1ClosedSearch());
	}
}
