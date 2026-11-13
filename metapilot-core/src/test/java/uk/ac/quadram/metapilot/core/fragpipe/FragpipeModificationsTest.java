package uk.ac.quadram.metapilot.core.fragpipe;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;

class FragpipeModificationsTest {

	@TempDir
	Path tempDir;

	private static List<String> rows(String table) {
		return Arrays.asList(table.split("; "));
	}

	@Test
	@DisplayName("Default fixed modification puts carbamidomethyl on C and leaves every other slot at zero")
	void defaultFixedTable() {
		List<String> rows = rows(FragpipeModifications.fixedTable(FragpipeModifications.DEFAULT_FIXED));
		assertEquals(30, rows.size());
		assertTrue(rows.contains("57.02146,C (cysteine),true,-1"));
		assertTrue(rows.contains("0.0,M (methionine),true,-1"));
		assertEquals(29, rows.stream().filter(r -> r.startsWith("0.0,")).count());
	}

	@Test
	@DisplayName("Default variable modifications are oxidation on M and protein N-terminal acetylation")
	void defaultVariableTable() {
		List<String> rows = rows(FragpipeModifications.variableTable(FragpipeModifications.DEFAULT_VARIABLE));
		assertEquals(FragpipeModifications.VARIABLE_SLOTS, rows.size());
		assertEquals("15.9949,M,true,3", rows.get(0));
		assertEquals("42.0106,[^,true,1", rows.get(1));
		assertEquals("0.0,site_3,false,1", rows.get(2));
		assertEquals("0.0,site_16,false,1", rows.get(15));
	}

	@Test
	@DisplayName("Terminus names map to MSFragger site codes and MSFragger's own syntax keeps its case")
	void siteCodes() {
		List<String> rows = rows(FragpipeModifications.variableTable(List.of("Carbamyl:peptide n-term:43.0058",
				"Amidated:protein c-term:-0.984", "Pyro-glu:nQnC:-17.0265", "Phospho:STY:79.96633")));
		assertEquals("43.0058,n^,true,1", rows.get(0));
		assertEquals("-0.984,]^,true,1", rows.get(1));
		assertEquals("-17.0265,nQnC,true,3", rows.get(2));
		assertEquals("79.96633,STY,true,3", rows.get(3));

		List<String> fixed = rows(FragpipeModifications.fixedTable(List.of("peptide n-term:229.1629", "k:229.1629")));
		assertTrue(fixed.contains("229.1629,N-Term Peptide,true,-1"));
		assertTrue(fixed.contains("229.1629,K (lysine),true,-1"));
	}

	@Test
	@DisplayName("Entries FragPipe cannot hold are rejected instead of silently dropped")
	void invalidEntriesAreRejected() {
		assertThrows(IllegalArgumentException.class, () -> FragpipeModifications.fixedTable(List.of("CM:57.02")));
		assertThrows(IllegalArgumentException.class, () -> FragpipeModifications.fixedTable(List.of("C=57.02")));
		assertThrows(IllegalArgumentException.class,
				() -> FragpipeModifications.variableTable(List.of("Oxidation:M:abc")));
		List<String> tooMany = new ArrayList<>();
		for (int i = 0; i <= FragpipeModifications.VARIABLE_SLOTS; i++) {
			tooMany.add("Mod" + i + ":M:1.0");
		}
		assertThrows(IllegalArgumentException.class, () -> FragpipeModifications.variableTable(tooMany));
	}

	@Test
	@DisplayName("apply replaces both tables of a workflow and leaves other keys alone")
	void applyReplacesWorkflowTables() throws Exception {
		Path workflow = tempDir.resolve("test.workflow");
		Files.writeString(workflow, "msfragger.table.fix-mods=0.0,C-Term Peptide,true,-1\n"
				+ "msfragger.table.var-mods=79.96633,STY,true,3\n" + "msfragger.search_enzyme_name_1=trypsin\n");
		FragpipeWorkflowFile file = FragpipeWorkflowFile.parse(workflow);

		FragpipeModifications.apply(file, FragpipeModifications.DEFAULT_FIXED, FragpipeModifications.DEFAULT_VARIABLE);

		assertEquals(FragpipeModifications.fixedTable(FragpipeModifications.DEFAULT_FIXED),
				file.get("msfragger.table.fix-mods"));
		assertTrue(file.get("msfragger.table.var-mods").startsWith("15.9949,M,true,3; 42.0106,[^,true,1; "));
		assertEquals("trypsin", file.get("msfragger.search_enzyme_name_1"));
	}
}
