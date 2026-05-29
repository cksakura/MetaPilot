package uk.ac.quadram.metapilot.task.par;

import uk.ac.quadram.metapilot.task.dia.DeepComprehensivePaths;

import java.io.File;
import java.util.Arrays;
import java.util.Comparator;
import java.util.List;
import java.util.Locale;

/**
 * Discovers DIA/DDA raw inputs from a folder for headless execution.
 */
public final class HpcInputDiscovery {

    private static final String[] DEFAULT_EXTENSIONS = { ".raw.dia", ".raw.d", ".dia", ".raw", ".mzml", ".d" };

    private HpcInputDiscovery() {
    }

    public static MetaData discoverMetaData(String inputDir) {
        if (inputDir == null || inputDir.isBlank()) {
            return null;
        }
        File dir = new File(inputDir.trim());
        if (!dir.isDirectory()) {
            throw new IllegalArgumentException("input_dir is not a directory: " + dir.getAbsolutePath());
        }
        File[] files = dir.listFiles(f -> f != null && f.isFile() && matchesExtension(f.getName()));
        if (files == null || files.length == 0) {
            throw new IllegalArgumentException("No supported raw files found under input_dir: " + dir.getAbsolutePath());
        }
        Arrays.sort(files, Comparator.comparing(f -> f.getName().toLowerCase(Locale.ROOT)));
        String[] rawFiles = new String[files.length];
        String[] expNames = new String[files.length];
        int[] fractions = new int[files.length];
        int[] replicates = new int[files.length];
        for (int i = 0; i < files.length; i++) {
            rawFiles[i] = files[i].getAbsolutePath();
            expNames[i] = DeepComprehensivePaths.runIdFromRawFilePath(rawFiles[i]);
            fractions[i] = 1;
            replicates[i] = 1;
        }
        return new MetaData(rawFiles, expNames, fractions, replicates);
    }

    private static boolean matchesExtension(String fileName) {
        String lower = fileName.toLowerCase(Locale.ROOT);
        for (String ext : DEFAULT_EXTENSIONS) {
            if (lower.endsWith(ext)) {
                return true;
            }
        }
        return false;
    }
}
