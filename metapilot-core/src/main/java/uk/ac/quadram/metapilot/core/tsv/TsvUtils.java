package uk.ac.quadram.metapilot.core.tsv;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.HashSet;

public final class TsvUtils {

    private TsvUtils() {
    }

    /**
     * Writes a TSV containing rows from {@code file1} whose first-column value does NOT appear
     * in the first column of {@code file2}.
     *
     * <p>By default, the first line of each input is treated as a header and is not considered a key.
     * The output will include {@code file1}'s header if present.</p>
     */
    public static void writeRowsFromFile1NotInFile2FirstColumn(File file1, File file2, File outFile) throws IOException {
        writeRowsFromFile1NotInFile2FirstColumn(file1, true, file2, true, outFile);
    }

    public static void writeRowsFromFile1NotInFile2FirstColumn(
            File file1,
            boolean file1HasHeader,
            File file2,
            boolean file2HasHeader,
            File outFile) throws IOException {

        if (file1 == null || !file1.isFile()) {
            throw new IOException("file1 not found: " + (file1 == null ? "null" : file1.getAbsolutePath()));
        }
        if (file2 == null || !file2.isFile()) {
            throw new IOException("file2 not found: " + (file2 == null ? "null" : file2.getAbsolutePath()));
        }
        if (outFile == null) {
            throw new IOException("outFile is null");
        }
        File parent = outFile.getParentFile();
        if (parent != null && !parent.exists() && !parent.mkdirs()) {
            throw new IOException("Cannot create output folder: " + parent.getAbsolutePath());
        }

        HashSet<String> file2Keys = new HashSet<>();
        try (BufferedReader r2 = new BufferedReader(new FileReader(file2))) {
            String line;
            if (file2HasHeader) {
                r2.readLine(); // header
            }
            while ((line = r2.readLine()) != null) {
                if (line.isEmpty()) continue;
                String key = firstColumn(line);
                if (!key.isEmpty()) {
                    file2Keys.add(key);
                }
            }
        }

        try (BufferedReader r1 = new BufferedReader(new FileReader(file1));
                PrintWriter w = new PrintWriter(outFile)) {
            String line;
            if (file1HasHeader) {
                String header = r1.readLine();
                if (header == null) {
                    return;
                }
                w.println(header);
            }
            while ((line = r1.readLine()) != null) {
                if (line.isEmpty()) continue;
                String key = firstColumn(line);
                if (key.isEmpty()) continue;
                if (!file2Keys.contains(key)) {
                    w.println(line);
                }
            }
        }
    }
    
    /**
     * Writes a TSV containing rows from {@code file1} whose first-column value does NOT appear
     * in the first column of {@code file2}.
     *
     * <p>By default, the first line of each input is treated as a header and is not considered a key.
     * The output will include {@code file1}'s header if present.</p>
     */
    public static void writeRowsFromFile1NotInFile2FirstColumn(File file1, File file2, File outFile, HashSet<String> set) throws IOException {
        writeRowsFromFile1NotInFile2FirstColumn(file1, true, file2, true, outFile, set);
    }

    public static void writeRowsFromFile1NotInFile2FirstColumn(
            File file1,
            boolean file1HasHeader,
            File file2,
            boolean file2HasHeader,
            File outFile, HashSet<String> set) throws IOException {

        if (file1 == null || !file1.isFile()) {
            throw new IOException("file1 not found: " + (file1 == null ? "null" : file1.getAbsolutePath()));
        }
        if (file2 == null || !file2.isFile()) {
            throw new IOException("file2 not found: " + (file2 == null ? "null" : file2.getAbsolutePath()));
        }
        if (outFile == null) {
            throw new IOException("outFile is null");
        }
        File parent = outFile.getParentFile();
        if (parent != null && !parent.exists() && !parent.mkdirs()) {
            throw new IOException("Cannot create output folder: " + parent.getAbsolutePath());
        }

        HashSet<String> file2Keys = new HashSet<>();
        try (BufferedReader r2 = new BufferedReader(new FileReader(file2))) {
            String line;
            if (file2HasHeader) {
                r2.readLine(); // header
            }
            while ((line = r2.readLine()) != null) {
                if (line.isEmpty()) continue;
                String key = firstColumn(line).replace('L', 'I');
                if(set.contains(key)) {
                	if (!key.isEmpty()) {
                        file2Keys.add(key);
                    }
                }
            }
        }
        System.out.println("file2Keys\t"+file2Keys.size());

        try (BufferedReader r1 = new BufferedReader(new FileReader(file1));
                PrintWriter w = new PrintWriter(outFile)) {
            String line;
            if (file1HasHeader) {
                String header = r1.readLine();
                if (header == null) {
                    return;
                }
                w.println(header);
            }
            while ((line = r1.readLine()) != null) {
                if (line.isEmpty()) continue;
                String key = firstColumn(line).replace('L', 'I');
                if (key.isEmpty()) continue;
                if(set.contains(key)) {
                	if (!file2Keys.contains(key)) {
                        w.println(line);
                    }
                }
            }
        }
    }

    private static String firstColumn(String tsvLine) {
        if (tsvLine == null) return "";
        int tab = tsvLine.indexOf('\t');
        String first = tab >= 0 ? tsvLine.substring(0, tab) : tsvLine;
        return first == null ? "" : first.trim();
    }
    
	public static void main(String[] args) throws IOException {
		File parentFile = new File("Z:\\PRIDE_dataset\\PXD051104_metaExpertPro\\DIA_combined\\MGYG000002293");
		HashSet<String> set = new HashSet<String>();
		try (BufferedReader r2 = new BufferedReader(
				new FileReader(new File(parentFile, "MGYG000002293.pr_matrix.tsv")))) {
			String line = r2.readLine();
			while ((line = r2.readLine()) != null) {
				if (line.isEmpty())
					continue;
				String[] cs = line.split("\t");
				set.add(cs[6].replace('L', 'I'));
			}
		}
		System.out.println(set.size());

		TsvUtils.writeRowsFromFile1NotInFile2FirstColumn(
				new File(parentFile, "MGYG000002293_unipept_6.5.0_UniProtKB_2025.04.tsv"), true,
				new File(parentFile, "Common_unipept_6.5.0_UniProtKB_2025.04.tsv"), true,
				new File(parentFile, "MetaPilot_unipept_6.5.0_UniProtKB_2025.04.tsv"), set);
	}
}

