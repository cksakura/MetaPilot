package uk.ac.quadram.metapilot.core.fragpipe;

public final class FragpipeWorkflowPaths {

    private FragpipeWorkflowPaths() {
    }

    public static String formatDatabasePath(String path) {
        if (path == null) return "";

        String normalized = path.trim().replace("\\\\", "\\").replaceFirst("^([A-Za-z])\\\\:", "$1:");
        if (normalized.isEmpty()) return "";

        if (isWindowsPath(normalized)) {
            return normalized.replace('/', '\\').replace("\\", "\\\\");
        }

        return normalized.replace('\\', '/');
    }

    private static boolean isWindowsPath(String path) {
        return path.matches("^[A-Za-z]:[\\\\/].*") || path.startsWith("\\");
    }
}
