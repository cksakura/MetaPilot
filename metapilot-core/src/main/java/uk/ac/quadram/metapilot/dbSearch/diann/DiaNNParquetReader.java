package uk.ac.quadram.metapilot.dbSearch.diann;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Set;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;

import org.apache.commons.io.FileUtils;
import org.apache.hadoop.conf.Configuration;
import org.apache.hadoop.fs.Path;
import org.apache.parquet.ParquetReadOptions;
import org.apache.parquet.column.page.PageReadStore;
import org.apache.parquet.example.data.Group;
import org.apache.parquet.example.data.simple.convert.GroupRecordConverter;
import org.apache.parquet.hadoop.ParquetFileReader;
import org.apache.parquet.hadoop.metadata.ParquetMetadata;
import org.apache.parquet.hadoop.util.HadoopInputFile;
import org.apache.parquet.io.ColumnIOFactory;
import org.apache.parquet.io.MessageColumnIO;
import org.apache.parquet.io.RecordReader;
import org.apache.parquet.io.InputFile;
import org.apache.parquet.schema.MessageType;
import org.apache.parquet.schema.PrimitiveType;
import org.apache.parquet.schema.Type;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import uk.ac.quadram.metapilot.core.aminoacid.AminoAcidProperty;

public class DiaNNParquetReader {

	private static final Logger LOGGER = LogManager.getLogger(DiaNNParquetReader.class);

	private DiaNNPrecursor[] diaNNPrecursors;

	static {
		// Avoid "HADOOP_HOME and hadoop.home.dir are unset" / winutils.exe warnings on Windows
		// when using Parquet (Hadoop libs) without a real Hadoop install
		if (System.getProperty("hadoop.home.dir") == null && System.getenv("HADOOP_HOME") == null) {
			String dir = System.getProperty("java.io.tmpdir");
			if (dir != null && !dir.isEmpty()) {
				System.setProperty("hadoop.home.dir", dir);
			}
		}
	}

	/** Silence verbose Parquet and Hadoop DEBUG logging (works with Logback or Log4j2). */
	private static void silenceParquetAndHadoopLogging() {
		try {
			// Try Logback (e.g. when slf4j binds to logback-classic)
			Object factory = Class.forName("org.slf4j.LoggerFactory").getMethod("getILoggerFactory").invoke(null);
			if (factory != null && "ch.qos.logback.classic.LoggerContext".equals(factory.getClass().getName())) {
				Object levelWarn = Class.forName("ch.qos.logback.classic.Level").getField("WARN").get(null);
				Class<?> levelClass = Class.forName("ch.qos.logback.classic.Level");
				for (String name : new String[] { "org.apache.parquet", "org.apache.hadoop" }) {
					Object log = factory.getClass().getMethod("getLogger", String.class).invoke(factory, name);
					if (log != null) log.getClass().getMethod("setLevel", levelClass).invoke(log, levelWarn);
				}
				return;
			}
		} catch (Throwable ignored) {
		}
		try {
			// Try Log4j2
			Class<?> logManager = Class.forName("org.apache.logging.log4j.LogManager");
			Class<?> level = Class.forName("org.apache.logging.log4j.Level");
			Object warn = level.getField("WARN").get(null);
			for (String name : new String[] { "org.apache.parquet", "org.apache.hadoop" }) {
				Object log = logManager.getMethod("getLogger", String.class).invoke(null, name);
				if (log != null) log.getClass().getMethod("setLevel", level).invoke(log, warn);
			}
		} catch (Throwable ignored) {
		}
	}
	
	public static HashMap<String, HashSet<String>> getPepProMap(String parquet) throws IOException {
		LOGGER.debug("Reading parquet file {}", parquet);
		
		HashMap<String, HashSet<String>> map = new HashMap<String, HashSet<String>>();
		Path path = new Path(parquet);

		HadoopInputFile inputFile = HadoopInputFile.fromPath(path, new Configuration());
		ParquetReadOptions readOptions = ParquetReadOptions.builder().build();
		ParquetFileReader fileReader = ParquetFileReader.open(inputFile, readOptions);

		// Get file schema
		ParquetMetadata metadata = fileReader.getFooter();
		MessageType schema = metadata.getFileMetaData().getSchema();

		// Loop through row groups in the Parquet file
		PageReadStore pages;
		while ((pages = fileReader.readNextRowGroup()) != null) {
			long rowCount = pages.getRowCount();
			LOGGER.debug("Parquet row group has {} rows", rowCount);

			MessageColumnIO columnIO = new ColumnIOFactory().getColumnIO(schema);
			GroupRecordConverter recordConverter = new GroupRecordConverter(schema);
			RecordReader<Group> recordReader = columnIO.getRecordReader(pages, recordConverter);

			for (int i = 0; i < rowCount; i++) {
				Group group = recordReader.read();
				String sequence = group.getBinary("Stripped.Sequence", 0).toStringUsingUTF8();
				String[] proteins = group.getBinary("Protein.Ids", 0).toStringUsingUTF8().split(";");
				HashSet<String> set;
				if (map.containsKey(sequence)) {
					set = map.get(sequence);
				} else {
					set = new HashSet<String>();
					map.put(sequence, set);
				}
				for (int j = 0; j < proteins.length; j++) {
					set.add(proteins[j]);
				}
			}
		}

		fileReader.close();

		return map;
	}
	
	public static HashMap<String, HashSet<String>> getProPepMap(String parquet) throws IOException {
		LOGGER.debug("Reading parquet file {}", parquet);

		HashMap<String, HashSet<String>> map = new HashMap<String, HashSet<String>>();
		Path path = new Path(parquet);

		HadoopInputFile inputFile = HadoopInputFile.fromPath(path, new Configuration());
		ParquetReadOptions readOptions = ParquetReadOptions.builder().build();
		ParquetFileReader fileReader = ParquetFileReader.open(inputFile, readOptions);

		// Get file schema
		ParquetMetadata metadata = fileReader.getFooter();
		MessageType schema = metadata.getFileMetaData().getSchema();

		// Loop through row groups in the Parquet file
		PageReadStore pages;
		while ((pages = fileReader.readNextRowGroup()) != null) {
			long rowCount = pages.getRowCount();
			LOGGER.debug("Parquet row group has {} rows", rowCount);

			MessageColumnIO columnIO = new ColumnIOFactory().getColumnIO(schema);
			GroupRecordConverter recordConverter = new GroupRecordConverter(schema);
			RecordReader<Group> recordReader = columnIO.getRecordReader(pages, recordConverter);

			for (int i = 0; i < rowCount; i++) {
				Group group = recordReader.read();
				String sequence = group.getBinary("Stripped.Sequence", 0).toStringUsingUTF8();
				String[] proteins = group.getBinary("Protein.Ids", 0).toStringUsingUTF8().split(";");
				for (int j = 0; j < proteins.length; j++) {
					HashSet<String> set;
					if (map.containsKey(proteins[j])) {
						set = map.get(proteins[j]);
					} else {
						set = new HashSet<String>();
						map.put(proteins[j], set);
					}
					set.add(sequence);
				}
			}
		}

		fileReader.close();
		return map;
	}

	public static HashSet<String> getPeptideSet(String parquet) throws IOException {
		LOGGER.debug("Reading parquet file {}", parquet);
		HashSet<String> set = new HashSet<String>();
		Path path = new Path(parquet);

		HadoopInputFile inputFile = HadoopInputFile.fromPath(path, new Configuration());
		ParquetReadOptions readOptions = ParquetReadOptions.builder().build();
		ParquetFileReader fileReader = ParquetFileReader.open(inputFile, readOptions);

		// Get file schema
		ParquetMetadata metadata = fileReader.getFooter();
		MessageType schema = metadata.getFileMetaData().getSchema();

		// Loop through row groups in the Parquet file
		PageReadStore pages;
		while ((pages = fileReader.readNextRowGroup()) != null) {
			long rowCount = pages.getRowCount();
			LOGGER.debug("Parquet row group has {} rows", rowCount);

			MessageColumnIO columnIO = new ColumnIOFactory().getColumnIO(schema);
			GroupRecordConverter recordConverter = new GroupRecordConverter(schema);
			RecordReader<Group> recordReader = columnIO.getRecordReader(pages, recordConverter);

			for (int i = 0; i < rowCount; i++) {
				Group group = recordReader.read();
				if ((int) group.getLong("Decoy", 0) != 1) {
					String sequence = group.getBinary("Stripped.Sequence", 0).toStringUsingUTF8();
					set.add(sequence);
				}
			}
		}

		fileReader.close();

		return set;
	}

	public static HashMap<String, Double> getPeptideMap(String parquet) throws IOException {
		HashMap<String, Double> map = new HashMap<String, Double>();
		LOGGER.debug("Reading parquet file {}", parquet);
		Path path = new Path(parquet);

		HadoopInputFile inputFile = HadoopInputFile.fromPath(path, new Configuration());
		ParquetReadOptions readOptions = ParquetReadOptions.builder().build();
		ParquetFileReader fileReader = ParquetFileReader.open(inputFile, readOptions);

		// Get file schema
		ParquetMetadata metadata = fileReader.getFooter();
		MessageType schema = metadata.getFileMetaData().getSchema();

		// Loop through row groups in the Parquet file
		PageReadStore pages;
		while ((pages = fileReader.readNextRowGroup()) != null) {
			long rowCount = pages.getRowCount();
			LOGGER.debug("Parquet row group has {} rows", rowCount);

			MessageColumnIO columnIO = new ColumnIOFactory().getColumnIO(schema);
			GroupRecordConverter recordConverter = new GroupRecordConverter(schema);
			RecordReader<Group> recordReader = columnIO.getRecordReader(pages, recordConverter);

			for (int i = 0; i < rowCount; i++) {
				Group group = recordReader.read();
				if ((int) group.getLong("Decoy", 0) != 1) {
					String sequence = group.getBinary("Stripped.Sequence", 0).toStringUsingUTF8();
					double preNormal = group.getFloat("Precursor.Normalised", 0);
					if (map.containsKey(sequence)) {
						map.put(sequence, map.get(sequence) + preNormal);
					} else {
						map.put(sequence, preNormal);
					}
				}
			}
		}

		fileReader.close();

		return map;
	}

	public DiaNNParquetReader(String parquet, String cut) {
		silenceParquetAndHadoopLogging();
		ArrayList<DiaNNPrecursor> list = new ArrayList<DiaNNPrecursor>();
		HashSet<Character> cutSet = new HashSet<Character>();
		HashSet<Character> noCutSet = new HashSet<Character>();
		String[] cuts = cut.split(",");
		for (int i = 0; i < cuts.length; i++) {
			if (cuts[i].startsWith("!*")) {
				noCutSet.add(cuts[i].charAt(2));
			} else {
				cutSet.add(cuts[i].charAt(0));
			}
		}

		try {
			Path path = new Path(parquet);
			LOGGER.debug("Reading parquet file {}", parquet);
			HadoopInputFile inputFile = HadoopInputFile.fromPath(path, new Configuration());
			ParquetReadOptions readOptions = ParquetReadOptions.builder().build();
			ParquetFileReader fileReader = ParquetFileReader.open(inputFile, readOptions);

			// Get file schema
			ParquetMetadata metadata = fileReader.getFooter();
			MessageType schema = metadata.getFileMetaData().getSchema();
			boolean hasMs1Area = false;
			boolean hasCScore = false;
			for (Type f : schema.getFields()) {
				if ("Ms1.Area".equals(f.getName())) {
					hasMs1Area = true;
				} else if ("CScore".equals(f.getName())) {
					hasCScore = true;
				}
			}
			// Print all column names
			StringBuilder colNames = new StringBuilder("Parquet columns: ");
			for (Type f : schema.getFields()) {
				if (colNames.length() > 16) colNames.append(", ");
				colNames.append(f.getName());
			}
			LOGGER.debug(colNames.toString());

			// Loop through row groups in the Parquet file
			PageReadStore pages;
			while ((pages = fileReader.readNextRowGroup()) != null) {
				long rowCount = pages.getRowCount();
//				System.out.println("Row group has " + rowCount + " rows.");

				MessageColumnIO columnIO = new ColumnIOFactory().getColumnIO(schema);
				GroupRecordConverter recordConverter = new GroupRecordConverter(schema);
				RecordReader<Group> recordReader = columnIO.getRecordReader(pages, recordConverter);

				for (int i = 0; i < rowCount; i++) {
					Group group = recordReader.read();
					String run = group.getBinary("Run", 0).toStringUsingUTF8();
					if (run.lastIndexOf(".") > 0) {
						run = run.substring(0, run.lastIndexOf("."));
					}
					String[] proteins = group.getBinary("Protein.Ids", 0).toStringUsingUTF8().split(";");
					String modSeq = group.getBinary("Modified.Sequence", 0).toStringUsingUTF8();
					String sequence = group.getBinary("Stripped.Sequence", 0).toStringUsingUTF8();
					int missCleavage = 0;
					for (int j = 0; j < sequence.length() - 1; j++) {
						char aa = sequence.charAt(j);
						char bb = sequence.charAt(j + 1);
						if (cutSet.contains(aa) && !noCutSet.contains(bb)) {
							missCleavage++;
						}
					}
					double precursorMz = group.getFloat("Precursor.Mz", 0);
					int charge = (int) group.getLong("Precursor.Charge", 0);
					// Precursor.Mz from DIA-NN is the charge-dependent m/z, not the neutral
					// peptide mass. Convert to neutral monoisotopic mass: M = (m/z - proton) * z
					// (same formula used in PixelList.java / PeakSegment.java elsewhere in this
					// codebase). Without this conversion, the "Mass" field silently held m/z.
					double mass = charge > 0 ? (precursorMz - AminoAcidProperty.PROTON_W) * charge : precursorMz;
					double preNormal = group.getFloat("Precursor.Normalised", 0);
					double preQuan = group.getFloat("Precursor.Quantity", 0);
					double ms1Normal = group.getFloat("Ms1.Normalised", 0);
					double ms1Area = hasMs1Area ? group.getFloat("Ms1.Area", 0) : 0;
					double qvalue = group.getFloat("Q.Value", 0);
					double PEP = group.getFloat("PEP", 0);
					double globalQ = group.getFloat("Global.Q.Value", 0);
					double pgQvalue = group.getFloat("PG.Q.Value", 0);
					double globalPGQvalue = group.getFloat("Global.PG.Q.Value", 0);
					boolean isDecoy = ((int) group.getLong("Decoy", 0)) == 1;
					DiaNNPrecursor diaNNPrecursor = new DiaNNPrecursor(run, proteins, modSeq, sequence, charge,
							missCleavage, mass, ms1Normal, qvalue, PEP, globalQ, pgQvalue, globalPGQvalue,
							preNormal, preQuan, isDecoy);
					if (hasMs1Area) {
						diaNNPrecursor.setMs1Area(ms1Area);
					}
					if (hasCScore) {
						double cscore = group.getFloat("CScore", 0);
						diaNNPrecursor.setCscore(cscore);
					}

					list.add(diaNNPrecursor);
				}
			}
			fileReader.close();
		} catch (IOException e) {
			LOGGER.error("Failed to read DIA-NN parquet file: " + parquet, e);
		}

		this.diaNNPrecursors = list.toArray(new DiaNNPrecursor[list.size()]);
	}

	public DiaNNPrecursor[] getDiaNNPrecursors() {
		return diaNNPrecursors;
	}

	public DiaNNPrecursor[] getDiaNNPrecursors(double globalQValue) {
		ArrayList<DiaNNPrecursor> list = new ArrayList<DiaNNPrecursor>();
		for (DiaNNPrecursor pp : diaNNPrecursors) {
			if (pp.getGlobalQvalue() < globalQValue) {
				list.add(pp);
			}
		}
		return list.toArray(DiaNNPrecursor[]::new);
	}
	
	public static ArrayList<String> getParquetColumnNames(String filePath) throws IOException {
		ArrayList<String> columnNames = new ArrayList<>();
		Path path = new Path(filePath);
		Configuration conf = new Configuration();

    // 1. Adapt the Hadoop Path to the generic InputFile interface
		InputFile inputFile = HadoopInputFile.fromPath(path, conf);

    // 2. Open the reader (this automatically reads the footer)
    // Using try-with-resources to ensure it closes automatically
		try (ParquetFileReader reader = ParquetFileReader.open(inputFile)) {
        
        // 3. Extract the schema (MessageType)
			MessageType schema = reader.getFileMetaData().getSchema();

        // 4. Iterate and collect names
			for (Type field : schema.getFields()) {
				columnNames.add(field.getName());
			}
		}

		return columnNames;
    }

	/**
	 * Export parquet contents to a TSV file for direct inspection.
	 * Writes one header line (column names) then one data row per parquet row.
	 * Cell values are tab-separated; tabs and newlines in values are replaced by spaces.
	 *
	 * @param parquetPath path to the parquet file
	 * @param tsvPath     path for the output TSV file
	 * @throws IOException if read/write fails
	 */
	public static void exportToTsv(String parquetPath, String tsvPath) throws IOException {
		exportToTsv(parquetPath, new File(tsvPath));
	}

	/**
	 * Export parquet contents to a TSV file for direct inspection.
	 *
	 * @param parquetPath path to the parquet file
	 * @param tsvFile     output TSV file
	 * @throws IOException if read/write fails
	 */
	public static void exportToTsv(String parquetPath, File tsvFile) throws IOException {
		silenceParquetAndHadoopLogging();
		Path path = new Path(parquetPath);
		HadoopInputFile inputFile = HadoopInputFile.fromPath(path, new Configuration());
		ParquetReadOptions readOptions = ParquetReadOptions.builder().build();
		ParquetFileReader fileReader = ParquetFileReader.open(inputFile, readOptions);

		ParquetMetadata metadata = fileReader.getFooter();
		MessageType schema = metadata.getFileMetaData().getSchema();
		int numFields = schema.getFieldCount();
		String[] columnNames = new String[numFields];
		for (int i = 0; i < numFields; i++) {
			columnNames[i] = schema.getFieldName(i);
		}

		try (PrintWriter out = new PrintWriter(tsvFile, StandardCharsets.UTF_8.name())) {
			// Header
			out.println(String.join("\t", columnNames));

			PageReadStore pages;
			while ((pages = fileReader.readNextRowGroup()) != null) {
				long rowCount = pages.getRowCount();
				MessageColumnIO columnIO = new ColumnIOFactory().getColumnIO(schema);
				GroupRecordConverter recordConverter = new GroupRecordConverter(schema);
				RecordReader<Group> recordReader = columnIO.getRecordReader(pages, recordConverter);

				for (long r = 0; r < rowCount; r++) {
					Group group = recordReader.read();
					String[] cells = new String[numFields];
					for (int i = 0; i < numFields; i++) {
						cells[i] = cellToString(schema, group, i);
					}
					out.println(String.join("\t", cells));
				}
			}
		} finally {
			fileReader.close();
		}
	}

	/** Format one cell value for TSV (escape tabs and newlines). */
	private static String cellToString(MessageType schema, Group group, int fieldIndex) {
		String raw;
		try {
			if (group.getFieldRepetitionCount(fieldIndex) == 0) {
				return "";
			}
			Type type = schema.getType(fieldIndex);
			if (type.isPrimitive()) {
				PrimitiveType pt = type.asPrimitiveType();
				switch (pt.getPrimitiveTypeName()) {
					case BINARY:
					case FIXED_LEN_BYTE_ARRAY:
						raw = group.getBinary(fieldIndex, 0).toStringUsingUTF8();
						break;
					case INT32:
						raw = String.valueOf(group.getInteger(fieldIndex, 0));
						break;
					case INT64:
						raw = String.valueOf(group.getLong(fieldIndex, 0));
						break;
					case FLOAT:
						raw = String.valueOf(group.getFloat(fieldIndex, 0));
						break;
					case DOUBLE:
						raw = String.valueOf(group.getDouble(fieldIndex, 0));
						break;
					case BOOLEAN:
						raw = String.valueOf(group.getBoolean(fieldIndex, 0));
						break;
					default:
						raw = String.valueOf(group.getValueToString(fieldIndex, 0));
				}
			} else {
				raw = String.valueOf(group.getValueToString(fieldIndex, 0));
			}
		} catch (Exception e) {
			raw = "";
		}
		// Keep TSV valid: replace tab and newline in values
		if (raw == null) return "";
		return raw.replace("\t", " ").replace("\n", " ").replace("\r", " ");
	}

	@SuppressWarnings("unused")
	private static void compare(File folder) {
		HashSet<String> prSet = new HashSet<String>();
		try (BufferedReader reader = new BufferedReader(
				new FileReader(new File(folder, folder.getName() + ".pr_matrix.tsv")))) {
			String line = reader.readLine();
			while ((line = reader.readLine()) != null) {
				String[] cs = line.split("\t");
				prSet.add(cs[6]);
			}
			reader.close();
		} catch (IOException e) {

		}

		File parquetFile = new File(folder, folder.getName() + ".parquet");
		DiaNNParquetReader reader = new DiaNNParquetReader(parquetFile.getAbsolutePath(), "K*R*");
		DiaNNPrecursor[] pps = reader.getDiaNNPrecursors();
		HashSet<String> ppSet = new HashSet<String>();
		HashSet<String> trSet = new HashSet<String>();
		for (DiaNNPrecursor pp : pps) {
			ppSet.add(pp.getSeqString());
			if (!pp.getIsDecoy()) {
				trSet.add(pp.getSeqString());
			}
		}

		HashSet<String> allSet = new HashSet<String>();
		allSet.addAll(trSet);
		allSet.addAll(prSet);
		int common = trSet.size() + prSet.size() - allSet.size();
		LOGGER.debug("{}\t{}\t{}\t{}", ppSet.size(), prSet.size(), common, trSet.size());
	}

	@SuppressWarnings("unused")
	private static void matchGenome(String dbFile, String parquetFile) {

		DiaNNParquetReader ppreader = new DiaNNParquetReader(parquetFile, "K*R*");
		DiaNNPrecursor[] pps = ppreader.getDiaNNPrecursors();
		HashSet<String> allPepSet = new HashSet<String>();
		for (DiaNNPrecursor pp : pps) {
			if (pp.getIsDecoy()) {
				allPepSet.add(pp.getSeqString());
			}
		}
		LOGGER.debug("decoy\t{}\t{}", pps.length, allPepSet.size());

		ConcurrentHashMap<String, HashMap<String, HashSet<String>>> genomePepMap = new ConcurrentHashMap<String, HashMap<String, HashSet<String>>>();
		HashSet<String> matchedPepSet = new HashSet<String>();
		ExecutorService executor = Executors.newFixedThreadPool(12);
		File currentFile = new File(dbFile);
		File predictPepFile = new File(currentFile, "predicted");
		File[] files = predictPepFile.listFiles();
		for (File file : files) {
			executor.execute(() -> {
				if (file.getName().endsWith(".predicted.tsv")) {
					try (BufferedReader reader = new BufferedReader(new FileReader(file))) {
						HashMap<String, HashSet<String>> proPepMap = new HashMap<String, HashSet<String>>();
						String pline = reader.readLine();
						while ((pline = reader.readLine()) != null) {
							String[] cs = pline.split("\t");
							if (allPepSet.contains(cs[1])) {
								matchedPepSet.add(cs[1]);
								if (proPepMap.containsKey(cs[0])) {
									proPepMap.get(cs[0]).add(cs[1]);
								} else {
									HashSet<String> pepSet = new HashSet<String>();
									pepSet.add(cs[1]);
									proPepMap.put(cs[0], pepSet);
								}
							}
						}
						if (proPepMap.size() > 0) {
							String preFileName = file.getName();
							String genome = preFileName.substring(0, preFileName.length() - ".predicted.tsv".length());
							genomePepMap.put(genome, proPepMap);
						}
					} catch (IOException e) {

					}
				}
			});
		}

		executor.shutdown();

		try {
			if (executor.awaitTermination(1, TimeUnit.HOURS)) {

				allPepSet.clear();
				allPepSet.addAll(matchedPepSet);

			} else {
			}
		} catch (InterruptedException e) {
		}

		LOGGER.debug("decoy\t{}\t{}", genomePepMap.size(), matchedPepSet.size());
	}
	
	private static void changeName(String in) {
		File[] files = (new File(in)).listFiles();
		for (int i = 0; i < files.length; i++) {
			String name = files[i].getName();
			if (name.endsWith(".parquet.parquet")) {
				String newName = name.substring(0, name.length() - ".parquet".length());
				files[i].renameTo(new File(in, newName));
			}
		}
	}
	
	public static void main(String[] args) throws IOException {
//		DiaNNParquetReader.changeName("Z:\\PRIDE_dataset\\PXD051104_metaExpertPro\\MetaPilot_DIA_MPW\\MetaPilot\\marker_deep\\marker");
//		DiaNNParquetReader
//				.getProPepMap("Z:\\PRIDE_dataset\\PXD051104_metaExpertPro\\MetaPilot_DIA_MPW\\MetaPilot\\marker_deep\\marker"
//						+ "\\per_file_1.parquet");
//		System.out.println("peptide count\t"+set.size());
//		DiaNNParquetReader.matchGenome("Z:\\Kai\\Database\\human-gut\\v2.0.2", 
//				"E:\\20250423\\20250422_MetaLAB_DIA_search\\selfModel\\mag_result\\combined.parquet");
		
		DiaNNParquetReader reader = new DiaNNParquetReader("Z:\\PRIDE_dataset\\IPX0003851000\\IPX0003851001\\DIA\\MetaPilot\\seed_deep\\seed\\round1"
				+ "\\round1-first-pass.parquet", "KR");
		DiaNNPrecursor[] pps =  reader.getDiaNNPrecursors(0.01);
		int[] count = new int[2];
		HashSet<String> set = new HashSet<String>();
		for(DiaNNPrecursor pp: pps) {
			if(pp.getIsDecoy()) {
				count[0]++;
			}else {
				count[1]++;
				set.add(pp.getSeqString());
			}
		}
		System.out.println(count[0]+"\t"+count[1]+"\t"+set.size());
		
		reader = new DiaNNParquetReader("Z:\\PRIDE_dataset\\IPX0003851000\\IPX0003851001\\DIA\\MetaPilot\\seed_deep\\seed\\round1"
				+ "\\round1.parquet", "KR");
		pps =  reader.getDiaNNPrecursors(0.01);
		count = new int[2];
		for(DiaNNPrecursor pp: pps) {
			if(pp.getIsDecoy()) {
				count[0]++;
			}else {
				count[1]++;
				set.add(pp.getSeqString());
			}
		}
		System.out.println(count[0]+"\t"+count[1]+"\t"+set.size());
		
//		System.out.println(DiaNNParquetReader.getPepProMap(
//				"C:\\Users\\daz25hak\\MetaProAI\\database\\MAG\\human-gut\\v2.0.2\\libraries\\Pekka"
//						+ "\\Pekka.parquet"));
		/*
		File[] files = (new File("Z:\\PRIDE_dataset\\IPX0003851000\\IPX0003851004\\DIA\\MetaPilot\\marker_deep\\marker")).listFiles();
		for (int i = 0; i < files.length; i++) {
			String name = files[i].getName();
			if(name.equals("marker")) {
				continue;
			}
			File parquetFile = new File(files[i], name + ".parquet");
			if (parquetFile.exists()) {
				File tsvFile = new File(files[i], name + ".parquet.tsv");
				if (!tsvFile.exists()) {
					DiaNNParquetReader.exportToTsv(parquetFile.getAbsolutePath(), tsvFile.getAbsolutePath());					
				}
				if (tsvFile.exists()) {
					FileUtils.copyFile(tsvFile,
							new File("Z:\\PRIDE_dataset\\IPX0003851000\\DATA_analysis\\data", tsvFile.getName()));
					Set<String> peptides = new HashSet<>();
					try {
						DiaNNParquetReader reader = new DiaNNParquetReader(parquetFile.getAbsolutePath(), "KR");
						DiaNNPrecursor[] precursors = reader.getDiaNNPrecursors();
						if (precursors != null) {
							for (DiaNNPrecursor p : precursors) {
								if (p == null || p.getIsDecoy()) continue;
								String seq = p.getSeqString();
								if (seq != null && !seq.isEmpty()) peptides.add(seq);
							}
						}
					} catch (Exception e) {

					}
					System.out.println(parquetFile.getName() + "\t" + peptides.size());
				}
			}
		}
		*/
	}

}
