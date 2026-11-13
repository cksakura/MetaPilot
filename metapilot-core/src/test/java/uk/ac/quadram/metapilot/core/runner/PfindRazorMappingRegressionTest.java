package uk.ac.quadram.metapilot.core.runner;

import static org.junit.jupiter.api.Assertions.*;
import org.junit.jupiter.api.Test;
import uk.ac.quadram.metapilot.quant.flashLFQ.FlashLfqQuanPeptide;
import uk.ac.quadram.metapilot.task.io.pep.MetaPeptide;

class PfindRazorMappingRegressionTest {
    @Test void modifiedPeptideCannotOverwriteUnmodifiedRazor() {
        MetaPeptide unmodified = new FlashLfqQuanPeptide("LAIADILVK", "LAIADILVK",
                new String[]{"MGYG000000119_00087"}, new double[]{1}, new int[]{1}, false);
        MetaPeptide modified = new FlashLfqQuanPeptide("LAIADILVK", "LAIADILVK(Lys->Arg)",
                new String[]{"MGYG000000123_02935"}, new double[]{1}, new int[]{1}, false);
        var map = PfindMagPostAnalysisTsvExporter.buildRazorProteinMap(
                new MetaPeptide[]{unmodified, modified},
                new String[]{"MGYG000000119_00087", "MGYG000000123_02935"});
        assertEquals("MGYG000000119_00087", map.get(unmodified.getModSeq()));
        assertEquals("MGYG000000123_02935", map.get(modified.getModSeq()));
        // A razor present in inference remains represented by its exact peptide.
        assertEquals(2, new java.util.HashSet<>(map.values()).size());
    }
}
