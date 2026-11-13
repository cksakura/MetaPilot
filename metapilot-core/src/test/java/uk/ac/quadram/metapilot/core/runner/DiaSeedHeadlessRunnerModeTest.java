package uk.ac.quadram.metapilot.core.runner;

import static org.junit.jupiter.api.Assertions.assertDoesNotThrow;
import static org.junit.jupiter.api.Assertions.assertThrows;

import java.io.IOException;

import org.junit.jupiter.api.Test;

import uk.ac.quadram.metapilot.task.dia.DiaRefinementMode;

/** Seed modes the headless DIA Seed runner accepts: Direct, Core / Focus, Comprehensive and Deep together. */
class DiaSeedHeadlessRunnerModeTest {

	@Test
	void acceptsTheFourSeedModesOfTheGui() {
		assertDoesNotThrow(() -> DiaSeedHeadlessRunner.requireSupportedMode(DiaRefinementMode.DIRECT));
		assertDoesNotThrow(() -> DiaSeedHeadlessRunner.requireSupportedMode(DiaRefinementMode.CORE));
		assertDoesNotThrow(() -> DiaSeedHeadlessRunner.requireSupportedMode(DiaRefinementMode.COMPREHENSIVE));
		assertDoesNotThrow(() -> DiaSeedHeadlessRunner.requireSupportedMode(DiaRefinementMode.DEEP_TOGETHER));
	}

	@Test
	void rejectsDeepPerFileAndMissingMode() {
		assertThrows(IOException.class, () -> DiaSeedHeadlessRunner.requireSupportedMode(DiaRefinementMode.DEEP_PER_FILE));
		assertThrows(IOException.class, () -> DiaSeedHeadlessRunner.requireSupportedMode(null));
	}

	@Test
	void persistedSeedModeValuesMapToTheRunnerModes() {
		assertDoesNotThrow(() -> DiaSeedHeadlessRunner.requireSupportedMode(DiaRefinementMode.fromSeedMode("DIRECT")));
		assertDoesNotThrow(() -> DiaSeedHeadlessRunner.requireSupportedMode(DiaRefinementMode.fromSeedMode("FOCUS")));
		assertDoesNotThrow(() -> DiaSeedHeadlessRunner.requireSupportedMode(DiaRefinementMode.fromSeedMode("COMPREHENSIVE")));
		assertDoesNotThrow(() -> DiaSeedHeadlessRunner.requireSupportedMode(DiaRefinementMode.fromSeedMode("DEEP")));
	}
}
