package uk.ac.quadram.metapilot.core.exec;

import java.util.ArrayList;
import java.util.List;

/**
 * Minimal command-line tokenizer (space-separated, supports quotes).
 * Intended for splitting parameter strings like: {@code --out "a b.txt" --threads 4}.
 */
public final class CommandLine {

    private CommandLine() {
    }

    public static List<String> split(String cmd) {
        List<String> out = new ArrayList<>();
        if (cmd == null) return out;

        String s = cmd.trim();
        if (s.isEmpty()) return out;

        StringBuilder cur = new StringBuilder();
        boolean inSingle = false;
        boolean inDouble = false;

        for (int i = 0; i < s.length(); i++) {
            char c = s.charAt(i);

            if (c == '\'' && !inDouble) {
                inSingle = !inSingle;
                continue;
            }
            if (c == '"' && !inSingle) {
                inDouble = !inDouble;
                continue;
            }

            if (!inSingle && !inDouble && Character.isWhitespace(c)) {
                if (cur.length() > 0) {
                    out.add(cur.toString());
                    cur.setLength(0);
                }
                continue;
            }

            cur.append(c);
        }

        if (cur.length() > 0) {
            out.add(cur.toString());
        }

        return out;
    }
}

