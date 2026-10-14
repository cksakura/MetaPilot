package uk.ac.quadram.metapilot.core.runner;

import static org.junit.jupiter.api.Assertions.assertArrayEquals;
import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertSame;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.io.File;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;

import uk.ac.quadram.metapilot.core.workflow.ProgressCallback;
import uk.ac.quadram.metapilot.spectra.io.MSConvertor;
import uk.ac.quadram.metapilot.task.dda.mag.par.DDAMagOpenSearchParameters;
import uk.ac.quadram.metapilot.task.par.MetaData;

class RawSpectraConversionTest {

	private static final long T0 = 1_700_000_000_000L;
	private static final String COMPLETE = "<indexedmzML>...</indexedmzML>\n";

	@TempDir
	Path tmp;

	@Test
	void mzmlInputIsLeftAsLoaded() throws IOException {
		File a = write(new File(tmp.toFile(), "a.mzML"), COMPLETE, T0);
		DDAMagOpenSearchParameters params = params(a);
		List<String[]> before = params.getExperimentData();

		assertSame(before, RawSpectraConversion.prepare(params, new ProgressCallback.NoOp()));
		assertEquals(a.getPath(), params.getExperimentData().get(0)[0]);
		assertEquals(a.getPath(), params.getMetadata().getRawFiles()[0]);
		assertFalse(new File(new File(params.getResult()), "spectra").exists());
	}

	@Test
	void mixedInputIsRefused() throws IOException {
		DDAMagOpenSearchParameters params = params(write(new File(tmp.toFile(), "a.raw"), "x", T0),
				write(new File(tmp.toFile(), "b.mzML"), COMPLETE, T0));

		IOException e = assertThrows(IOException.class,
				() -> RawSpectraConversion.prepare(params, new ProgressCallback.NoOp()));
		assertTrue(e.getMessage().contains("b.mzML"), e.getMessage());
	}

	@Test
	void completeConversionInSpectraIsReusedWithoutMsconvert() throws IOException {
		File raw = write(new File(tmp.toFile(), "S1.RAW"), "x", T0);
		DDAMagOpenSearchParameters params = params(raw);
		File mzml = write(new File(params.getSpectraFile(), "S1.mzML"), COMPLETE, T0 + 1000);

		List<String[]> rows = RawSpectraConversion.prepare(params, new ProgressCallback.NoOp());

		assertEquals(mzml.getAbsolutePath(), rows.get(0)[0]);
		assertEquals("S1", rows.get(0)[1]);
		assertSame(rows, params.getExperimentData());
		assertArrayEquals(new String[] { mzml.getAbsolutePath() }, params.getMetadata().getRawFiles());
	}

	@Test
	void truncatedConversionIsNotReused() throws IOException {
		File raw = write(new File(tmp.toFile(), "S1.raw"), "x", T0);
		DDAMagOpenSearchParameters params = params(raw);
		write(new File(params.getSpectraFile(), "S1.mzML"), "<indexedmzML><mzML>", T0 + 1000);

		IOException e = assertThrows(IOException.class,
				() -> RawSpectraConversion.prepare(params, new ProgressCallback.NoOp()));
		assertTrue(e.getMessage().contains("msconvert is not configured"), e.getMessage());
	}

	@Test
	void conversionOlderThanItsRawIsNotReused() throws IOException {
		File raw = write(new File(tmp.toFile(), "S1.raw"), "x", T0 + 1000);
		DDAMagOpenSearchParameters params = params(raw);
		write(new File(params.getSpectraFile(), "S1.mzML"), COMPLETE, T0);

		IOException e = assertThrows(IOException.class,
				() -> RawSpectraConversion.prepare(params, new ProgressCallback.NoOp()));
		assertTrue(e.getMessage().contains("msconvert is not configured"), e.getMessage());
		assertEquals(raw.getPath(), params.getExperimentData().get(0)[0]);
	}

	@Test
	void sameSampleNameTwiceIsRefused() throws IOException {
		File a = write(new File(new File(tmp.toFile(), "x"), "S1.raw"), "x", T0);
		File b = write(new File(new File(tmp.toFile(), "y"), "S1.raw"), "x", T0);

		IOException e = assertThrows(IOException.class,
				() -> RawSpectraConversion.prepare(params(a, b), new ProgressCallback.NoOp()));
		assertTrue(e.getMessage().contains("S1"), e.getMessage());
	}

	@Test
	void recognisesAnUnknownInstrument() {
		assertTrue(RawSpectraConversion.isVersionProblem(Arrays.asList("processing file: x.raw",
				"[Reader_Thermo::fillInMetadata] unable to parse instrument model; make sure you are using the "
						+ "latest version of ProteoWizard")));
		assertFalse(RawSpectraConversion.isVersionProblem(Arrays.asList("writing output file: x.mzML")));
	}

	@Test
	void commandMatchesRaw2mzML() {
		File raw = new File("C:\\data dir\\S1.raw");
		File out = new File("C:\\out dir");
		assertEquals(Arrays.asList("msconvert.exe", raw.getAbsolutePath(), "--mzML", "-o", out.getAbsolutePath(),
				"--zlib", "--filter", "peakPicking true 1-"), MSConvertor.raw2mzMLCommand("msconvert.exe", raw, out));
	}

	private DDAMagOpenSearchParameters params(File... spectra) {
		DDAMagOpenSearchParameters params = new DDAMagOpenSearchParameters();
		File result = new File(tmp.toFile(), "result");
		result.mkdirs();
		params.setResult(result.getAbsolutePath());
		List<String[]> rows = new ArrayList<>();
		String[] files = new String[spectra.length];
		String[] names = new String[spectra.length];
		for (int i = 0; i < spectra.length; i++) {
			files[i] = spectra[i].getPath();
			names[i] = spectra[i].getName().substring(0, spectra[i].getName().lastIndexOf('.'));
			rows.add(new String[] { files[i], names[i], "1" });
		}
		params.setExperimentData(rows);
		params.setMetadata(new MetaData(files, names, new int[spectra.length], new int[spectra.length]));
		return params;
	}

	private static File write(File file, String content, long modified) throws IOException {
		file.getParentFile().mkdirs();
		Files.write(file.toPath(), content.getBytes(StandardCharsets.UTF_8));
		file.setLastModified(modified);
		return file;
	}
}
