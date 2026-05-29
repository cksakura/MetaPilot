package uk.ac.quadram.metapilot.core.runner;

import uk.ac.quadram.metapilot.core.workflow.ProgressCallback;
import uk.ac.quadram.metapilot.task.io.pro.MetaProtein;
import uk.ac.quadram.metapilot.task.mag.MagDbItem;
import uk.ac.quadram.metapilot.task.par.MetaConstants;

import java.io.File;
import java.io.IOException;

/**
 * Headless MAG functional annotation export.
 */
public final class DiaMagFunctionalAnnotationExporter {

    private DiaMagFunctionalAnnotationExporter() {
    }

    public static void export(
            File resultDir,
            MagDbItem magDb,
            File funcDbFile,
            String microDbCatalogueId,
            String quantitativeMode,
            MetaProtein[] proteins,
            String[] intensityTitles,
            ProgressCallback cb
    ) throws Exception {
        if (cb == null) cb = new ProgressCallback.NoOp();
        if (resultDir == null) throw new IOException("resultDir is null");
        if (magDb == null) throw new IOException("magDb is null");
        if (funcDbFile == null || !funcDbFile.isFile()) {
            throw new IOException("Function definition database not found: "
                    + (funcDbFile == null ? "null" : funcDbFile.getAbsolutePath()));
        }
        if (proteins == null || proteins.length == 0) {
            cb.onMessage("DIA functional annotation: no proteins to annotate");
            return;
        }
        if (microDbCatalogueId == null || microDbCatalogueId.isBlank()) {
            throw new IOException("MAG catalogue id is required for functional annotation");
        }

        cb.onMessage("DIA functional annotation: annotating proteins");
        MagFunctionalAnnotationService.Request request = new MagFunctionalAnnotationService.Request();
        request.resultDir = resultDir;
        request.magDb = magDb;
        request.microDbCatalogueId = microDbCatalogueId;
        request.functionDefinitionDb = funcDbFile;
        request.analysisType = MetaConstants.pFindMag;
        request.quantitativeMode = quantitativeMode;
        request.metaProteins = proteins;
        request.intensityTitles = intensityTitles;
        MagFunctionalAnnotationService.Result result = MagFunctionalAnnotationService.run(request);
        cb.onMessage("DIA functional annotation: wrote " + result.annotatedProteinsTsv.getName());
    }
}
