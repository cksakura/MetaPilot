package uk.ac.quadram.metapilot.task.io.qpx;

import java.io.File;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.time.OffsetDateTime;
import java.time.format.DateTimeFormatter;
import java.util.Collection;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

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
import org.apache.parquet.schema.MessageType;
import org.apache.parquet.schema.MessageTypeParser;

import uk.ac.quadram.metapilot.MetaPilotVersion;
import uk.ac.quadram.metapilot.core.tools.LocalParquetOutputFile;
import uk.ac.quadram.metapilot.task.io.pep.MetaPSM;

/**
 * Writes a QPX {@code psm.parquet} (spec 1.1, https://qpx.quantms.org/spec/psm/) from
 * any {@link MetaPSM} collection. It reads only MetaPSM's normalised accessors, so it
 * is the same code path for every DDA workflow (FragPipe, Sage, pFind, open search,
 * MaxQuant).
 * <p>
 * Column types follow the spec: int64 {@code psm_id}, int16 {@code charge} and
 * {@code missed_cleavages}, float32 m/z / RT / ppm, float64 PEP, lists as standard
 * 3-level Parquet LISTs. Optional columns MetaPilot cannot fill (cross-links, spectral
 * arrays, feature link) are present and null so the file has the full schema.
 * <p>
 * {@code psm_id} is supplied by the producer (QPX keeps a producer ID): a stable
 * 64-bit FNV-1a hash of peptidoform, charge, run, scan and engine, made unique within
 * the file. A PSM without the required values (sequence, charge, observed mass, run)
 * is skipped and counted, never written with invented values.
 */
public final class QpxPsmParquetWriter {

	private static final Logger LOGGER = LogManager.getLogger(QpxPsmParquetWriter.class);

	public static final String QPX_VERSION = "1.1";
	public static final String DEFAULT_FILE_NAME = "metapilot.psm.parquet";

	/** Proton mass (Da), for neutral mass -> m/z. */
	public static final double PROTON = 1.007276466812;

	static {
		// Same Windows workaround as TsvToParquetExporter / TypedParquetTableWriter.
		if (System.getProperty("hadoop.home.dir") == null && System.getenv("HADOOP_HOME") == null) {
			String dir = System.getProperty("java.io.tmpdir");
			if (dir != null && !dir.isEmpty()) {
				System.setProperty("hadoop.home.dir", dir);
			}
		}
	}

	static final String SCHEMA = "message psm {\n"
			+ "  required int64 psm_id;\n"
			+ "  required binary sequence (UTF8);\n"
			+ "  required binary peptidoform (UTF8);\n"
			+ "  optional group modifications (LIST) {\n"
			+ "    repeated group list {\n"
			+ "      optional group element {\n"
			+ "        required binary name (UTF8);\n"
			+ "        optional binary accession (UTF8);\n"
			+ "        optional group positions (LIST) {\n"
			+ "          repeated group list {\n"
			+ "            optional group element {\n"
			+ "              required int32 position;\n"
			+ "              optional binary amino_acid (UTF8);\n"
			+ "              optional group scores (LIST) {\n"
			+ "                repeated group list {\n"
			+ "                  optional group element {\n"
			+ "                    required binary score_name (UTF8);\n"
			+ "                    required float score_value;\n"
			+ "                    optional boolean higher_better;\n"
			+ "                  }\n"
			+ "                }\n"
			+ "              }\n"
			+ "            }\n"
			+ "          }\n"
			+ "        }\n"
			+ "      }\n"
			+ "    }\n"
			+ "  }\n"
			+ "  required int32 charge (INT_16);\n"
			+ "  optional double posterior_error_probability;\n"
			+ "  required boolean is_decoy;\n"
			+ "  optional float calculated_mz;\n"
			+ "  required float observed_mz;\n"
			+ "  optional float mass_error_ppm;\n"
			+ "  optional int32 missed_cleavages (INT_16);\n"
			+ "  optional float rt;\n"
			+ "  optional float predicted_rt;\n"
			+ "  required binary run_file_name (UTF8);\n"
			+ "  required group scan (LIST) {\n"
			+ "    repeated group list {\n"
			+ "      required int32 element;\n"
			+ "    }\n"
			+ "  }\n"
			+ "  optional group additional_scores (LIST) {\n"
			+ "    repeated group list {\n"
			+ "      optional group element {\n"
			+ "        required binary score_name (UTF8);\n"
			+ "        required float score_value;\n"
			+ "        optional boolean higher_better;\n"
			+ "      }\n"
			+ "    }\n"
			+ "  }\n"
			+ "  optional group cv_params (LIST) {\n"
			+ "    repeated group list {\n"
			+ "      optional group element {\n"
			+ "        required binary cv_name (UTF8);\n"
			+ "        optional binary cv_value (UTF8);\n"
			+ "      }\n"
			+ "    }\n"
			+ "  }\n"
			+ "  optional group protein_accessions (LIST) {\n"
			+ "    repeated group list {\n"
			+ "      optional binary element (UTF8);\n"
			+ "    }\n"
			+ "  }\n"
			+ "  optional group cross_links (LIST) {\n"
			+ "    repeated group list {\n"
			+ "      optional group element {\n"
			+ "        required binary xl_type (UTF8);\n"
			+ "        optional binary partner_sequence (UTF8);\n"
			+ "        optional binary partner_peptidoform (UTF8);\n"
			+ "        required int32 donor_position;\n"
			+ "        optional int32 acceptor_position;\n"
			+ "        required binary linker_name (UTF8);\n"
			+ "        optional binary linker_accession (UTF8);\n"
			+ "        required double linker_mass;\n"
			+ "      }\n"
			+ "    }\n"
			+ "  }\n"
			+ "  optional float ion_mobility;\n"
			+ "  optional group mz_array (LIST) { repeated group list { optional float element; } }\n"
			+ "  optional group intensity_array (LIST) { repeated group list { optional float element; } }\n"
			+ "  optional group charge_array (LIST) { repeated group list { optional int32 element; } }\n"
			+ "  optional group ion_type_array (LIST) { repeated group list { optional binary element (UTF8); } }\n"
			+ "  optional group ion_mobility_array (LIST) { repeated group list { optional float element; } }\n"
			+ "  optional int64 feature_id;\n"
			+ "}\n";

	/** What happened in one write. */
	public static final class Summary {
		public int written;
		public int skipped;
		public int decoys;
		public int withPep;
		public int modificationsUnparsed;
		public final Set<String> engines = new LinkedHashSet<>();
		public File file;

		@Override
		public String toString() {
			return written + " PSMs written" + (skipped > 0 ? ", " + skipped + " skipped (missing required values)" : "")
					+ ", " + decoys + " decoys, " + withPep + " with PEP, engines " + engines
					+ (modificationsUnparsed > 0 ? ", " + modificationsUnparsed + " with unparsed modifications" : "");
		}
	}

	private QpxPsmParquetWriter() {
	}

	/** Parsed schema (exposed for tests). */
	public static MessageType schema() {
		return MessageTypeParser.parseMessageType(SCHEMA);
	}

	public static Summary write(File out, Collection<? extends MetaPSM> psms) throws IOException {
		if (out == null) {
			throw new IOException("output file is null");
		}
		MessageType schema = schema();
		Configuration conf = new Configuration();
		GroupWriteSupport.setSchema(schema, conf);

		File parent = out.getParentFile();
		if (parent != null) {
			parent.mkdirs();
		}
		File tmp = new File(parent, out.getName() + ".tmp");
		if (tmp.exists()) {
			tmp.delete();
		}
		File crc = new File(parent, "." + tmp.getName() + ".crc");
		if (crc.exists()) {
			crc.delete();
		}

		Summary summary = new Summary();
		Set<Long> usedIds = new HashSet<>();
		Map<String, String> meta = new LinkedHashMap<>();
		meta.put("qpx_version", QPX_VERSION);
		meta.put("software_provider", "MetaPilot " + MetaPilotVersion.getCurrentVersion());
		meta.put("scan_format", "scan");
		meta.put("creator", "MetaPilot");
		meta.put("file_type", "psm_file");
		meta.put("creation_date", OffsetDateTime.now().format(DateTimeFormatter.ISO_OFFSET_DATE_TIME));
		meta.put("compression_format", "gzip");

		// Engines are only known after the pass, but the footer is written at close;
		// collect them first so the metadata can name them.
		if (psms != null) {
			for (MetaPSM psm : psms) {
				if (psm != null && psm.getSearchEngine() != null) {
					summary.engines.add(psm.getSearchEngine());
				}
			}
		}
		meta.put("metapilot:search_engines", String.join(",", summary.engines));
		meta.put("metapilot:is_decoy_meaning",
				"search-engine target/decoy; MetaPilot shuffled genome-FDR entries are targets here and are "
						+ "flagged by cv_param metapilot:genome_fdr_decoy");

		try (ParquetWriter<Group> writer = ExampleParquetWriter.builder(new LocalParquetOutputFile(tmp))
				.withConf(conf).withType(schema).withExtraMetaData(meta)
				.withCompressionCodec(CompressionCodecName.GZIP)
				.withWriteMode(ParquetFileWriter.Mode.OVERWRITE).build()) {
			SimpleGroupFactory factory = new SimpleGroupFactory(schema);
			if (psms != null) {
				for (MetaPSM psm : psms) {
					Group row = toRow(factory, psm, usedIds, summary);
					if (row == null) {
						summary.skipped++;
						continue;
					}
					writer.write(row);
					summary.written++;
				}
			}
		}

		if (out.exists() && !out.delete()) {
			throw new IOException("Cannot replace " + out.getAbsolutePath());
		}
		if (!tmp.renameTo(out)) {
			throw new IOException("Cannot move " + tmp.getAbsolutePath() + " to " + out.getAbsolutePath());
		}
		if (crc.exists()) {
			crc.delete();
		}
		summary.file = out;
		return summary;
	}

	static Group toRow(SimpleGroupFactory factory, MetaPSM psm, Set<Long> usedIds, Summary summary) {
		if (psm == null) {
			return null;
		}
		String sequence = psm.getBaseSequence();
		int charge = psm.getCharge();
		String run = psm.getSpectrumFileName();
		double observedMass = psm.getObservedNeutralMass();
		if (sequence == null || sequence.isEmpty() || charge <= 0 || run == null || run.isEmpty()
				|| Double.isNaN(observedMass) || observedMass <= 0) {
			return null;
		}

		List<QpxModification> mods = psm.getQpxModifications();
		if (mods == null) {
			summary.modificationsUnparsed++;
		}
		String peptidoform = psm.getProForma();
		double theoreticalMass = psm.getTheoreticalNeutralMass();
		double observedMz = (observedMass + charge * PROTON) / charge;
		double calculatedMz = Double.isNaN(theoreticalMass) || theoreticalMass <= 0 ? Double.NaN
				: (theoreticalMass + charge * PROTON) / charge;

		Group g = factory.newGroup();
		g.append("psm_id", uniqueId(psm, peptidoform, run, usedIds));
		g.append("sequence", sequence.toUpperCase(java.util.Locale.ROOT));
		g.append("peptidoform", peptidoform);

		if (mods != null) {
			Group modList = g.addGroup("modifications");
			for (QpxModification mod : mods) {
				Group el = modList.addGroup("list").addGroup("element");
				el.append("name", mod.getName());
				if (mod.getAccession() != null) {
					el.append("accession", mod.getAccession());
				}
				Group positions = el.addGroup("positions");
				for (QpxModification.Site site : mod.getSites()) {
					Group p = positions.addGroup("list").addGroup("element");
					p.append("position", site.position);
					if (site.aminoAcid != null) {
						p.append("amino_acid", site.aminoAcid);
					}
					p.addGroup("scores"); // defined, empty: no localization scores
				}
			}
		}

		g.append("charge", charge);
		double pep = psm.getPosteriorErrorProbability();
		if (!Double.isNaN(pep)) {
			g.append("posterior_error_probability", pep);
			summary.withPep++;
		}
		boolean decoy = psm.isDecoy();
		g.append("is_decoy", decoy);
		if (decoy) {
			summary.decoys++;
		}
		if (!Double.isNaN(calculatedMz)) {
			g.append("calculated_mz", (float) calculatedMz);
			g.append("mass_error_ppm", (float) (1e6 * (observedMz - calculatedMz) / calculatedMz));
		}
		g.append("observed_mz", (float) observedMz);
		int miss = psm.getMiss();
		if (miss >= 0) {
			g.append("missed_cleavages", miss);
		}
		appendFloatIfFinite(g, "rt", psm.getRtSeconds());
		appendFloatIfFinite(g, "predicted_rt", psm.getPredictedRtSeconds());
		g.append("run_file_name", run);
		g.addGroup("scan").addGroup("list").append("element", psm.getScan());

		List<QpxScore> scores = psm.getAdditionalScores();
		if (scores != null && !scores.isEmpty()) {
			Group list = g.addGroup("additional_scores");
			for (QpxScore s : scores) {
				Group el = list.addGroup("list").addGroup("element");
				el.append("score_name", s.name);
				el.append("score_value", (float) s.value);
				if (s.higherBetter != null) {
					el.append("higher_better", s.higherBetter.booleanValue());
				}
			}
		}

		Map<String, String> cv = new LinkedHashMap<>();
		if (psm.getSearchEngine() != null) {
			cv.put("metapilot:search_engine", psm.getSearchEngine());
		}
		if (psm.isGenomeFdrDecoy()) {
			cv.put("metapilot:genome_fdr_decoy", "true");
		}
		if (psm.getModSeq() != null) {
			// The key the PSM joins final_peptides.tsv "Sequence" on.
			cv.put("metapilot:peptide_key", psm.getModSeq());
		}
		if (mods == null) {
			cv.put("metapilot:modifications_unparsed", psm.getFullSequence());
		}
		Map<String, String> extra = psm.getCvParams();
		if (extra != null) {
			cv.putAll(extra);
		}
		if (!cv.isEmpty()) {
			Group list = g.addGroup("cv_params");
			for (Map.Entry<String, String> e : cv.entrySet()) {
				Group el = list.addGroup("list").addGroup("element");
				el.append("cv_name", e.getKey());
				if (e.getValue() != null) {
					el.append("cv_value", e.getValue());
				}
			}
		}

		String[] proteins = psm.getProteins();
		if (proteins != null && proteins.length > 0) {
			Group list = g.addGroup("protein_accessions");
			for (String p : proteins) {
				Group el = list.addGroup("list");
				if (p != null && !p.trim().isEmpty()) {
					el.append("element", p.trim());
				}
			}
		}
		appendFloatIfFinite(g, "ion_mobility", psm.getIonMobility());
		return g;
	}

	private static void appendFloatIfFinite(Group g, String field, double value) {
		if (!Double.isNaN(value) && !Double.isInfinite(value)) {
			g.append(field, (float) value);
		}
	}

	/** FNV-1a 64 over the identity composite plus engine, probed to be unique in the file. */
	static long uniqueId(MetaPSM psm, String peptidoform, String run, Set<Long> used) {
		String key = peptidoform + '|' + psm.getCharge() + '|' + run + '|' + psm.getScan() + '|'
				+ psm.getSearchEngine();
		long h = 0xcbf29ce484222325L;
		for (byte b : key.getBytes(StandardCharsets.UTF_8)) {
			h ^= (b & 0xff);
			h *= 0x100000001b3L;
		}
		while (!used.add(h)) {
			h++;
		}
		return h;
	}
}
