package uk.ac.quadram.metapilot.core.tools;

import org.apache.hadoop.conf.Configuration;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.parquet.example.data.Group;
import org.apache.parquet.example.data.simple.SimpleGroupFactory;
import org.apache.parquet.hadoop.ParquetFileWriter;
import org.apache.parquet.hadoop.ParquetWriter;
import org.apache.parquet.hadoop.example.ExampleParquetWriter;
import org.apache.parquet.hadoop.example.GroupWriteSupport;
import org.apache.parquet.hadoop.metadata.CompressionCodecName;
import org.apache.parquet.schema.LogicalTypeAnnotation;
import org.apache.parquet.schema.MessageType;
import org.apache.parquet.schema.PrimitiveType;
import org.apache.parquet.schema.Types;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

/**
 * Mirrors an already-written result TSV file into an additional Parquet file with the
 * exact same information: same column names (in order) and the same row values.
 *
 * This is intentionally a "same content, different container" export: it does not
 * attempt to infer numeric/boolean types or re-derive values, it just re-encodes the
 * text already produced for the TSV report as Parquet so downstream tools (DuckDB,
 * pandas/pyarrow, MetaAtlas, etc.) can read the same table without parsing TSV.
 *
 * Every column is written as an OPTIONAL UTF8 (string) column. Type-aware schemas
 * (matching the MetaAtlas Parquet Export Upgrade Guide, e.g. int32/float64 columns)
 * can be layered on top of this later without changing the TSV writers.
 *
 * Uses the parquet-hadoop / hadoop-client stack already present on the classpath
 * (see DiaNNParquetReader) - no new dependencies required.
 */
public final class TsvToParquetExporter {

    private static final Logger LOGGER = LogManager.getLogger(TsvToParquetExporter.class);

    static {
        // Avoid "HADOOP_HOME and hadoop.home.dir are unset" / winutils.exe warnings on Windows
        // when using Parquet (Hadoop libs) without a real Hadoop install.
        if (System.getProperty("hadoop.home.dir") == null && System.getenv("HADOOP_HOME") == null) {
            String dir = System.getProperty("java.io.tmpdir");
            if (dir != null && !dir.isEmpty()) {
                System.setProperty("hadoop.home.dir", dir);
            }
        }
    }

    private TsvToParquetExporter() {
    }

    /**
     * Converts {@code tsvFile} to a sibling Parquet file (same name, ".parquet" extension),
     * logging a warning and returning {@code null} instead of throwing if anything goes wrong.
     * Safe to call unconditionally after writing a TSV report - it will never fail the caller's
     * export.
     */
    public static File convertQuietly(File tsvFile) {
        if (tsvFile == null) {
            return null;
        }
        try {
            return convert(tsvFile);
        } catch (Exception e) {
            LOGGER.warn("Skipping Parquet export for {}: {}", tsvFile.getAbsolutePath(), e.getMessage());
            return null;
        }
    }

    /**
     * Converts {@code tsvFile} to a sibling Parquet file (same name, ".parquet" extension).
     */
    public static File convert(File tsvFile) throws IOException {
        File parquetFile = new File(tsvFile.getParentFile(), stripExtension(tsvFile.getName()) + ".parquet");
        convert(tsvFile, parquetFile);
        return parquetFile;
    }

    /**
     * Converts {@code tsvFile} to {@code parquetFile}, preserving header and row order.
     */
    public static void convert(File tsvFile, File parquetFile) throws IOException {
        if (tsvFile == null || !tsvFile.exists()) {
            throw new IOException("TSV file not found: " + (tsvFile == null ? "null" : tsvFile.getAbsolutePath()));
        }
        if (parquetFile == null) {
            throw new IOException("parquetFile is null");
        }

        String[] header;
        List<String[]> rows = new ArrayList<>();
        try (BufferedReader br = new BufferedReader(
                new InputStreamReader(new FileInputStream(tsvFile), StandardCharsets.UTF_8))) {
            String headerLine = br.readLine();
            header = (headerLine == null) ? new String[0] : headerLine.split("\t", -1);
            String line;
            while ((line = br.readLine()) != null) {
                if (line.isEmpty()) {
                    continue;
                }
                rows.add(line.split("\t", -1));
            }
        }

        if (header.length == 0) {
            LOGGER.warn("TSV file {} has no header/columns; skipping Parquet export", tsvFile.getAbsolutePath());
            return;
        }

        writeParquet(parquetFile, safeColumnNames(header), rows);
    }

    private static void writeParquet(File parquetFile, String[] columns, List<String[]> rows) throws IOException {
        // NOTE: intermediate .named(...) calls resolve through the inherited
        // BaseGroupBuilder<MessageType, GroupBuilder<MessageType>>, so they return
        // Types.GroupBuilder<MessageType> (not Types.MessageTypeBuilder) even though the
        // runtime object is the same builder. The loop variable must be typed accordingly;
        // only the final .named(...) call (below) yields the built MessageType.
        Types.GroupBuilder<MessageType> schemaBuilder = Types.buildMessage();
        for (String column : columns) {
            schemaBuilder = schemaBuilder
                    .optional(PrimitiveType.PrimitiveTypeName.BINARY)
                    .as(LogicalTypeAnnotation.stringType())
                    .named(column);
        }
        MessageType schema = schemaBuilder.named("tsv_row");

        Configuration conf = new Configuration();
        GroupWriteSupport.setSchema(schema, conf);

        File parent = parquetFile.getParentFile();
        if (parent != null) {
            parent.mkdirs();
        }
        if (parquetFile.exists()) {
            parquetFile.delete();
        }

        try (ParquetWriter<Group> writer = ExampleParquetWriter.builder(new LocalParquetOutputFile(parquetFile))
                .withConf(conf)
                .withType(schema)
                .withCompressionCodec(CompressionCodecName.GZIP)
                .withWriteMode(ParquetFileWriter.Mode.OVERWRITE)
                .build()) {
            SimpleGroupFactory factory = new SimpleGroupFactory(schema);
            for (String[] row : rows) {
                Group group = factory.newGroup();
                for (int i = 0; i < columns.length; i++) {
                    String value = (i < row.length) ? row[i] : null;
                    if (value != null && !value.isEmpty()) {
                        group.append(columns[i], value);
                    }
                }
                writer.write(group);
            }
        }
    }

    /**
     * Parquet requires non-null, non-empty, unique field names within a message type.
     * TSV headers are otherwise used as-is (including spaces/punctuation) so column
     * names match the TSV exactly whenever possible.
     */
    private static String[] safeColumnNames(String[] header) {
        String[] out = new String[header.length];
        Set<String> used = new HashSet<>();
        for (int i = 0; i < header.length; i++) {
            String name = (header[i] == null) ? "" : header[i].trim();
            if (name.isEmpty()) {
                name = "column_" + (i + 1);
            }
            String candidate = name;
            int suffix = 2;
            while (!used.add(candidate)) {
                candidate = name + "_" + suffix++;
            }
            out[i] = candidate;
        }
        return out;
    }

    private static String stripExtension(String fileName) {
        int dot = fileName.lastIndexOf('.');
        return (dot > 0) ? fileName.substring(0, dot) : fileName;
    }
}
