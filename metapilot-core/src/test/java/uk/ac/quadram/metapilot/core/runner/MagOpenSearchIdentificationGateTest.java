package uk.ac.quadram.metapilot.core.runner;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import org.junit.jupiter.api.AfterEach;
import org.junit.jupiter.api.Test;

import uk.ac.quadram.metapilot.core.workflow.ProgressCallback;

class MagOpenSearchIdentificationGateTest {

	@AfterEach
	void clearProperty() {
		System.clearProperty(MagOpenSearchHeadlessRunner.CONTINUE_WITH_SUCCESSFUL_PROPERTY);
	}

	private static List<String[]> samples(String... names) {
		List<String[]> out = new ArrayList<>();
		for (String name : names) {
			out.add(new String[] { name });
		}
		return out;
	}

	@Test
	void noSkippedSamplesPassesThrough() throws Exception {
		List<String[]> successful = samples("a.mzML");
		List<String[]> out = MagOpenSearchHeadlessRunner.applyIdentificationGate(successful,
				Collections.emptyList(), null, new ProgressCallback.NoOp());
		assertEquals(1, out.size());
	}

	@Test
	void allSkippedThrows() {
		IOException ex = assertThrows(IOException.class,
				() -> MagOpenSearchHeadlessRunner.applyIdentificationGate(Collections.emptyList(),
						List.of("empty_a", "empty_b"), null, new ProgressCallback.NoOp()));
		assertTrue(ex.getMessage().contains("No samples have stage-1"));
	}

	@Test
	void partialWithoutOptInThrows() {
		IOException ex = assertThrows(IOException.class,
				() -> MagOpenSearchHeadlessRunner.applyIdentificationGate(samples("ok.mzML"), List.of("empty_a"),
						null, new ProgressCallback.NoOp()));
		assertTrue(ex.getMessage().contains("--continue-with-successful-files"));
	}

	@Test
	void partialWithPropertyContinues() throws Exception {
		System.setProperty(MagOpenSearchHeadlessRunner.CONTINUE_WITH_SUCCESSFUL_PROPERTY, "true");
		List<String[]> out = MagOpenSearchHeadlessRunner.applyIdentificationGate(samples("ok.mzML", "ok2.mzML"),
				List.of("empty_a"), null, new ProgressCallback.NoOp());
		assertEquals(2, out.size());
	}

	@Test
	void partialWithHandlerCanContinue() throws Exception {
		List<String[]> out = MagOpenSearchHeadlessRunner.applyIdentificationGate(samples("ok.mzML"),
				List.of("empty_a"), (skipped, n) -> true, new ProgressCallback.NoOp());
		assertEquals(1, out.size());
	}

	@Test
	void partialWithHandlerCanStop() {
		IOException ex = assertThrows(IOException.class,
				() -> MagOpenSearchHeadlessRunner.applyIdentificationGate(samples("ok.mzML"), List.of("empty_a"),
						(skipped, n) -> false, new ProgressCallback.NoOp()));
		assertTrue(ex.getMessage().contains("Stopped after incomplete"));
	}
}
