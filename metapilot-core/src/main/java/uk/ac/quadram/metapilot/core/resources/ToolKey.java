package uk.ac.quadram.metapilot.core.resources;

import java.util.Locale;

/**
 * Stable keys used in meta_resources.params for CLI/headless mode.
 */
public enum ToolKey {
    JAVA,
    PYTHON,

    FRAGPIPE,
    FRAGPIPE_TOOLS,

    DIANN,
    DEEPDETECT,

    /** pFind.exe; its parent folder is the pFind bin directory (Windows only). */
    PFIND,
    /** FlashLFQ CMD.exe (Windows only for the pFind workflow). */
    FLASHLFQ,

    FUNCDB;

    public String key() {
        return name();
    }

    public static ToolKey fromStringOrNull(String s) {
        if (s == null) return null;
        String norm = s.trim().toUpperCase(Locale.ROOT);
        if (norm.isEmpty()) return null;
        try {
            return ToolKey.valueOf(norm);
        } catch (IllegalArgumentException e) {
            return null;
        }
    }
}

