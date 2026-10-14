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

import java.io.File;
import java.io.IOException;
import java.util.List;

/**
 * Writes a typed, "long" (one row per observation) Parquet table from in-memory rows -
 * real column types (double/int/boolean/string), one row per entity-per-sample, with a
 * sample column - rather than a wide, all-string mirror of an existing TSV.
 *
 * This is the building block for tables shaped to match the MetaAtlas Parquet Export
 * Upgrade Guide (e.g. taxon_abundance.parquet, function_abundance.parquet), which need
 * a sample column and typed abundance values to be queryable/joinable across projects
 * in DuckDB the way the guide's example queries assume.
 *
 * Uses the same parquet-hadoop / hadoop-client stack already on the classpath (see
 * DiaNNParquetReader) - no new dependencies required.
 */
public final class TypedParquetTableWriter {

    private static final Logger LOGGER = LogManager.getLogger(TypedParquetTableWriter.class);

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

    private TypedParquetTableWriter() {
    }

    public enum ColumnType {
        STRING, DOUBLE, INT, BOOLEAN
    }

    public static final class Column {
        public final String name;
        public final ColumnType type;

        public Column(String name, ColumnType type) {
            this.name = name;
            this.type = type;
        }
    }

    /**
     * Strips a leading "Intensity " prefix (used throughout MetaPilot's TSV headers, e.g.
     * "Intensity SampleA") to recover the underlying sample/experiment name. Note: this is
     * the raw experiment/run label as it appears in the search results - not a resolved,
     * SDRF-backed stable sample_accession per the MetaAtlas guide (section 12); that
     * resolution step is not yet implemented.
     */
    public static String[] sampleIdsFromIntensityTitles(String[] intensityTitles) {
        if (intensityTitles == null) {
            return new String[0];
        }
        String[] out = new String[intensityTitles.length];
        String prefix = "Intensity ";
        for (int i = 0; i < intensityTitles.length; i++) {
            String t = (intensityTitles[i] == null) ? "" : intensityTitles[i].trim();
            if (t.regionMatches(true, 0, prefix, 0, prefix.length())) {
                t = t.substring(prefix.length()).trim();
            }
            out[i] = t;
        }
        return out;
    }

    /**
     * Writes {@code rows} (each an Object[] aligned with {@code columns}) to {@code out} as
     * a Parquet file with a schema built from {@code columns}' declared types. Missing/null
     * cells are left unset (the columns are all OPTIONAL).
     */
    public static void write(File out, List<Column> columns, List<Object[]> rows) throws IOException {
        if (out == null) {
            throw new IOException("output file is null");
        }
        if (columns == null || columns.isEmpty()) {
            LOGGER.warn("No columns given; skipping Parquet export for {}", out.getAbsolutePath());
            return;
        }

        // NOTE: intermediate .named(...) calls resolve through the inherited
        // BaseGroupBuilder<MessageType, GroupBuilder<MessageType>>, so they return
        // Types.GroupBuilder<MessageType> (not Types.MessageTypeBuilder) even though the
        // runtime object is the same builder; only the final .named(...) call below yields
        // the built MessageType. The loop variable must be typed accordingly.
        Types.GroupBuilder<MessageType> schemaBuilder = Types.buildMessage();
        for (Column c : columns) {
            switch (c.type) {
                case STRING:
                    schemaBuilder = schemaBuilder
                            .optional(PrimitiveType.PrimitiveTypeName.BINARY)
                            .as(LogicalTypeAnnotation.stringType())
                            .named(c.name);
                    break;
                case DOUBLE:
                    schemaBuilder = schemaBuilder
                            .optional(PrimitiveType.PrimitiveTypeName.DOUBLE)
                            .named(c.name);
                    break;
                case INT:
                    schemaBuilder = schemaBuilder
                            .optional(PrimitiveType.PrimitiveTypeName.INT32)
                            .named(c.name);
                    break;
                case BOOLEAN:
                    schemaBuilder = schemaBuilder
                            .optional(PrimitiveType.PrimitiveTypeName.BOOLEAN)
                            .named(c.name);
                    break;
                default:
                    throw new IOException("Unsupported column type: " + c.type);
            }
        }
        MessageType schema = schemaBuilder.named("record");

        Configuration conf = new Configuration();
        GroupWriteSupport.setSchema(schema, conf);

        File parent = out.getParentFile();
        if (parent != null) {
            parent.mkdirs();
        }
        if (out.exists()) {
            out.delete();
        }

        try (ParquetWriter<Group> writer = ExampleParquetWriter.builder(new LocalParquetOutputFile(out))
                .withConf(conf)
                .withType(schema)
                .withCompressionCodec(CompressionCodecName.GZIP)
                .withWriteMode(ParquetFileWriter.Mode.OVERWRITE)
                .build()) {
            SimpleGroupFactory factory = new SimpleGroupFactory(schema);
            for (Object[] row : rows) {
                Group group = factory.newGroup();
                for (int i = 0; i < columns.size(); i++) {
                    Object value = (i < row.length) ? row[i] : null;
                    if (value == null) {
                        continue;
                    }
                    Column c = columns.get(i);
                    switch (c.type) {
                        case STRING:
                            String s = value.toString();
                            if (!s.isEmpty()) {
                                group.append(c.name, s);
                            }
                            break;
                        case DOUBLE:
                            group.append(c.name, ((Number) value).doubleValue());
                            break;
                        case INT:
                            group.append(c.name, ((Number) value).intValue());
                            break;
                        case BOOLEAN:
                            group.append(c.name, (Boolean) value);
                            break;
                        default:
                            break;
                    }
                }
                writer.write(group);
            }
        }
    }
}
