package uk.ac.quadram.metapilot.core.runner;

import java.io.File;
import java.io.IOException;

import org.apache.logging.log4j.Logger;

import uk.ac.quadram.metapilot.core.diann.DiaNNHeadlessRunner;
import uk.ac.quadram.metapilot.core.workflow.ProgressCallback;
import uk.ac.quadram.metapilot.dbSearch.diann.DiannParameter;
import uk.ac.quadram.metapilot.task.dia.DiaLibSearchPar;
import uk.ac.quadram.metapilot.task.dia.par.MetaParameterDiaMarker;

/**
 * Shared DIA-NN no-quant library search runner for DIA MAG workflows.
 */
public final class DiaNNNoQuantLibrarySearchService {

    private DiaNNNoQuantLibrarySearchService() {
    }

    public static void run(File diannExe, MetaParameterDiaMarker params, String[] diaFiles,
                           File libraryFile, File outParquet, int threads, boolean reanalyse,
                           boolean exportDecoy, ProgressCallback cb, File condaPrefix,
                           Logger logger, String taskName) throws Exception {
        if (libraryFile == null || !libraryFile.exists()) {
            throw new IOException("Library not found for DIA-NN no-quant search: "
                    + (libraryFile == null ? "null" : libraryFile.getAbsolutePath()));
        }
        if (outParquet == null) {
            throw new IOException("Output parquet is null for DIA-NN no-quant search");
        }
        if (diaFiles == null || diaFiles.length == 0) {
            throw new IOException("No DIA files supplied for DIA-NN no-quant search");
        }

        DiannParameter diannPar = DiaMarkerHeadlessRunner.buildDiannPar(params);
        DiaLibSearchPar libSearchPar = new DiaLibSearchPar();
        configureNoQuantSearch(diannPar, libSearchPar, diaFiles, threads, reanalyse, exportDecoy);

        String cmd = diannPar.generateLibSearchNoQuant(libSearchPar, outParquet.getAbsolutePath(),
                libraryFile.getAbsolutePath());
        cmd = cmd + DiaMarkerHeadlessRunner.buildFArgs(diaFiles);
        int exit = DiaNNHeadlessRunner.run(diannExe, cmd, outParquet.getParentFile(), cb, condaPrefix);
        if (exit != 0 && logger != null) {
            String name = taskName != null && !taskName.isBlank() ? taskName : "DIA-NN no-quant search";
            logger.warn("{}: DIA-NN no-quant search exit code {} (exportDecoy={}, out={})",
                    name, exit, exportDecoy, outParquet.getAbsolutePath());
        }
    }

    public static void configureNoQuantSearch(DiannParameter diannPar, DiaLibSearchPar libSearchPar,
                                              String[] diaFiles, int threads, boolean reanalyse,
                                              boolean exportDecoy) {
        if (diannPar == null) {
            throw new IllegalArgumentException("diannPar is null");
        }
        if (libSearchPar == null) {
            throw new IllegalArgumentException("libSearchPar is null");
        }
        int threadCount = Math.max(1, threads);
        libSearchPar.setMbr(false);
        libSearchPar.setDiaFiles(diaFiles);
        libSearchPar.setThreadCount(threadCount);
        libSearchPar.setqValue(0.01);
        diannPar.setThreads(threadCount);
        diannPar.setMatrices(false);
        diannPar.setRelaxed_prot_inf(true);
        diannPar.setReanalyse(reanalyse);
        diannPar.setExportDecoy(exportDecoy);
    }
}
