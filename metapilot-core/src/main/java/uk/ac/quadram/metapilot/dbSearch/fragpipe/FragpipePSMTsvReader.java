package uk.ac.quadram.metapilot.dbSearch.fragpipe;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Set;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import uk.ac.quadram.metapilot.core.prodb.DecoyShuffler;

/**
 * Reads a FragPipe psm.tsv file and returns a list of {@link FragPipePSM}.
 * Supports the standard FragPipe psm.tsv header (Spectrum, Spectrum File, Peptide,
 * Modified Peptide, Charge, Retention, Observed Mass, Calculated Peptide Mass,
 * Delta Mass, Hyperscore, Probability, Qvalue, Number of Missed Cleavages, Protein,
 * Assigned Modifications, Mapped Proteins, etc.). The delegate protein is in "Protein";
 * additional proteins are in "Mapped Proteins" (comma-separated). Both are combined for the
 * full protein list.
 */
public class FragpipePSMTsvReader {

	private static final Logger LOGGER = LogManager.getLogger(FragpipePSMTsvReader.class);

	/** Pattern for spectrum title like "file.123.456.789" (scan is first number group). */
	private static final Pattern TITLE_PATTERN = Pattern.compile("^.+\\.(\\d+)\\.(\\d+)\\.(\\d+)$");
	/** Fallback: "scan=123" style. */
	private static final Pattern SCAN_PATTERN = Pattern.compile("scan=(\\d+)", Pattern.CASE_INSENSITIVE);

	private final File psmFile;
	private String[] header;
	private int idxSpectrum;
	private int idxSpectrumFile;
	private int idxPeptide;
	private int idxModifiedPeptide;
	private int idxCharge;
	private int idxRetention;
	private int idxObservedMass;
	private int idxCalibratedObservedMass;
	private int idxCalculatedPeptideMass;
	private int idxDeltaMass;
	private int idxHyperscore;
	private int idxProbability;
	private int idxQvalue;
	private int idxIsDecoy;
	private int idxMissedCleavages;
	private int idxProtein;
	private int idxMappedProteins;
	private int idxAssignedModifications;

	public FragpipePSMTsvReader(File psmFile) {
		this.psmFile = psmFile;
	}

	public FragpipePSMTsvReader(String path) {
		this(new File(path));
	}

	/**
	 * Parse the header line and set column indices. Call once before reading rows.
	 */
	private void parseHeader(String headerLine) {
		header = headerLine.split("\t");
		idxSpectrum = -1;
		idxSpectrumFile = -1;
		idxPeptide = -1;
		idxModifiedPeptide = -1;
		idxCharge = -1;
		idxRetention = -1;
		idxObservedMass = -1;
		idxCalibratedObservedMass = -1;
		idxCalculatedPeptideMass = -1;
		idxDeltaMass = -1;
		idxHyperscore = -1;
		idxProbability = -1;
		idxQvalue = -1;
		idxIsDecoy = -1;
		idxMissedCleavages = -1;
		idxProtein = -1;
		idxMappedProteins = -1;
		idxAssignedModifications = -1;

		for (int i = 0; i < header.length; i++) {
			String col = header[i].trim();
			switch (col) {
				case "Spectrum":
					idxSpectrum = i;
					break;
				case "Spectrum File":
					idxSpectrumFile = i;
					break;
				case "Peptide":
					idxPeptide = i;
					break;
				case "Modified Peptide":
					idxModifiedPeptide = i;
					break;
				case "Charge":
					idxCharge = i;
					break;
				case "Retention":
					idxRetention = i;
					break;
				case "Observed Mass":
					idxObservedMass = i;
					break;
				case "Calibrated Observed Mass":
					idxCalibratedObservedMass = i;
					break;
				case "Calculated Peptide Mass":
					idxCalculatedPeptideMass = i;
					break;
				case "Delta Mass":
					idxDeltaMass = i;
					break;
				case "Hyperscore":
					idxHyperscore = i;
					break;
				case "Probability":
					idxProbability = i;
					break;
				case "Qvalue":
				case "Q-value":
					idxQvalue = i;
					break;
				case "Is Decoy":
				case "IsDecoy":
				case "Decoy":
					idxIsDecoy = i;
					break;
				case "Number of Missed Cleavages":
					idxMissedCleavages = i;
					break;
				case "Protein":
					idxProtein = i;
					break;
				case "Mapped Proteins":
					idxMappedProteins = i;
					break;
				case "Assigned Modifications":
					idxAssignedModifications = i;
					break;
				default:
					break;
			}
		}

		if (idxPeptide < 0 || idxProtein < 0) {
			LOGGER.warn("psm.tsv missing required columns (Peptide and/or Protein). Found columns: {}", (Object) header);
		}
	}

	/**
	 * Parse scan number from Spectrum column (e.g. "file.123.456.789" or "scan=123").
	 */
	public static int parseScanNumberOrDefault(String spectrum) {
		if (spectrum == null || spectrum.isEmpty())
			return 0;
		String s = spectrum.trim();
		Matcher m = TITLE_PATTERN.matcher(s);
		if (m.matches()) {
			try {
				long scan = Long.parseLong(m.group(1));
				if (scan > 0 && scan <= Integer.MAX_VALUE)
					return (int) scan;
			} catch (NumberFormatException ignored) {
				// fall through
			}
		}
		m = SCAN_PATTERN.matcher(s);
		if (m.find()) {
			try {
				long scan = Long.parseLong(m.group(1));
				if (scan > 0 && scan <= Integer.MAX_VALUE)
					return (int) scan;
			} catch (NumberFormatException ignored) {
				// fall through
			}
		}
		return 0;
	}

	private int maxRequiredIndex() {
		int max = 0;
		for (int i : new int[] { idxSpectrum, idxSpectrumFile, idxPeptide, idxModifiedPeptide, idxCharge,
				idxRetention, idxObservedMass, idxCalibratedObservedMass, idxCalculatedPeptideMass, idxDeltaMass,
				idxHyperscore, idxProbability, idxQvalue, idxIsDecoy, idxMissedCleavages, idxProtein, idxMappedProteins, idxAssignedModifications }) {
			if (i > max) max = i;
		}
		return max;
	}

	/**
	 * Read all PSMs from the file.
	 *
	 * @return list of FragPipePSM (never null)
	 * @throws IOException if the file cannot be read
	 */
	public List<FragPipePSM> read() throws IOException {
		List<FragPipePSM> list = new ArrayList<>();
		try (BufferedReader reader = new BufferedReader(new FileReader(psmFile))) {
			String line = reader.readLine();
			if (line == null) {
				LOGGER.warn("Empty psm.tsv: {}", psmFile.getAbsolutePath());
				return list;
			}
			parseHeader(line);
			int maxIdx = maxRequiredIndex();
			int lineNum = 1;
			while ((line = reader.readLine()) != null) {
				lineNum++;
				if (line.trim().isEmpty())
					continue;
				String[] cs = line.split("\t", -1);
				if (cs.length <= maxIdx) {
					LOGGER.trace("Skipping short line {} at line {}", lineNum, psmFile.getName());
					continue;
				}
				try {
					FragPipePSM psm = parseRow(cs);
					if (psm != null)
						list.add(psm);
				} catch (Exception e) {
					LOGGER.warn("Skipping malformed PSM at line {} in {}: {}", lineNum, psmFile.getName(), e.getMessage());
				}
			}
		}
		LOGGER.info("Read {} PSMs from {}", list.size(), psmFile.getAbsolutePath());
		return list;
	}

	/**
	 * Read all PSMs and return as array.
	 */
	public FragPipePSM[] getPsms() throws IOException {
		List<FragPipePSM> list = read();
		return list.toArray(new FragPipePSM[0]);
	}

	private String get(String[] row, int index) {
		if (index < 0 || index >= row.length)
			return "";
		String s = row[index];
		return s == null ? "" : s.trim();
	}

	private double getDouble(String[] row, int index, double defaultValue) {
		String s = get(row, index);
		if (s.isEmpty())
			return defaultValue;
		try {
			return Double.parseDouble(s);
		} catch (NumberFormatException e) {
			return defaultValue;
		}
	}

	private int getInt(String[] row, int index, int defaultValue) {
		String s = get(row, index);
		if (s.isEmpty())
			return defaultValue;
		try {
			return Integer.parseInt(s);
		} catch (NumberFormatException e) {
			return defaultValue;
		}
	}

	private FragPipePSM parseRow(String[] cs) {
		String spectrum = get(cs, idxSpectrum);
		String fileName = get(cs, idxSpectrumFile);
		String sequence = get(cs, idxPeptide);
		String modSeq = idxModifiedPeptide >= 0 ? get(cs, idxModifiedPeptide) : sequence;
		String mod = get(cs, idxAssignedModifications);
		int scan = parseScanNumberOrDefault(spectrum);
		int charge = getInt(cs, idxCharge, 0);
		double rt = getDouble(cs, idxRetention, 0.0);
		// Prefer Calibrated Observed Mass if present, else Observed Mass
		double expMass = idxCalibratedObservedMass >= 0 && !get(cs, idxCalibratedObservedMass).isEmpty()
				? getDouble(cs, idxCalibratedObservedMass, 0.0)
				: getDouble(cs, idxObservedMass, 0.0);
		double calMass = getDouble(cs, idxCalculatedPeptideMass, 0.0);
		double massDiff = getDouble(cs, idxDeltaMass, 0.0);
		int miss = getInt(cs, idxMissedCleavages, 0);
		double score = getDouble(cs, idxHyperscore, 0.0);
		double probability = getDouble(cs, idxProbability, 0.0);
		double qValue = getDouble(cs, idxQvalue, 1.0);
		boolean decoy = getBoolean(cs, idxIsDecoy, false);

		// Combine delegate protein (Protein) with Mapped Proteins for full protein list
		Set<String> proteinSet = new LinkedHashSet<>();
		String delegate = get(cs, idxProtein);
		if (!delegate.isEmpty()) {
			proteinSet.add(delegate.trim());
		}
		String mappedCol = idxMappedProteins >= 0 ? get(cs, idxMappedProteins) : "";
		if (!mappedCol.isEmpty()) {
			for (String p : mappedCol.split(",")) {
				String t = p.trim();
				if (!t.isEmpty()) {
					proteinSet.add(t);
				}
			}
		}
		String[] proteins = proteinSet.toArray(new String[0]);

		if (modSeq.isEmpty())
			modSeq = sequence;

		return new FragPipePSM(fileName, scan, spectrum, sequence, modSeq, mod, charge,
				expMass, calMass, massDiff, miss, rt, score, probability, qValue, proteins, decoy);
	}

	private boolean getBoolean(String[] row, int index, boolean defaultValue) {
		String s = get(row, index);
		if (s.isEmpty()) {
			return defaultValue;
		}
		return "true".equalsIgnoreCase(s) || "t".equalsIgnoreCase(s) || "yes".equalsIgnoreCase(s)
				|| "y".equalsIgnoreCase(s) || "1".equals(s);
	}
	
	private static void combineToFasta(String in, String host) {
		HashSet<String> set = new HashSet<String>();
		File[] files = (new File(in)).listFiles();
		for (int i = 0; i < files.length; i++) {
			String name = files[i].getName();
			name = name.replaceAll("-", "_");
			File subFile = new File(files[i], name + "_1_1");
			if (subFile.exists()) {
				File psmFile = new File(subFile, "psm.tsv");
				if (psmFile.exists()) {
					FragpipePSMTsvReader reader = new FragpipePSMTsvReader(psmFile);
					try {
						FragPipePSM[] psms = reader.getPsms();
						System.out.println(psms.length);
						for (int j = 0; j < psms.length; j++) {
							set.add(psms[j].getBaseSequence());
						}
					} catch (IOException e) {
						// TODO Auto-generated catch block
						e.printStackTrace();
					}
				}
			}
		}
		System.out.println("set\t" + set.size());

		File combinedFastaFile = new File(in, "peptide.fasta");
		try (PrintWriter writer = new PrintWriter(combinedFastaFile)) {
			int proid = 0;
			for (String seq : set) {
				String proName = "pro" + proid + " pro" + proid;
				writer.println(">" + proName);
				writer.println(seq);

				String shuffleProName = ">shuf_" + proName;
				String shuffle = DecoyShuffler.generateShuffledDecoy(seq, "", 42);
				writer.println(shuffleProName);
				writer.println(shuffle);

				String revProName = ">rev_" + proName;
				String revSeqString = (new StringBuilder(seq)).reverse().toString();

				writer.println(revProName);
				writer.println(revSeqString);

				String revShufProName = ">rev_shuf_" + proName;
				String revShufSeqString = (new StringBuilder(shuffle)).reverse().toString();

				writer.println(revShufProName);
				writer.println(revShufSeqString);

				proid++;
			}

			BufferedReader hostReader = new BufferedReader(new FileReader(host));
			String line = null;
			String hostProName = null;
			String refHostProName = null;
			StringBuilder seqSb = new StringBuilder();
			while ((line = hostReader.readLine()) != null) {
				if (line.startsWith(">")) {
					if (hostProName != null) {
						String seq = seqSb.toString();
						writer.println(hostProName);
						writer.println(seq);

						String cleavage = "KR";
						String shuffle = DecoyShuffler.generateShuffledDecoy(seq, cleavage, 42);
						writer.println(">shuf_" + hostProName.substring(1));
						writer.println(shuffle);

						writer.println(refHostProName);
						writer.println(new StringBuilder(seq).reverse().toString());

						writer.println(">rev_shuf_" + hostProName.substring(1));
						writer.println(new StringBuilder(shuffle).reverse().toString());
					}

					hostProName = line;
					refHostProName = ">rev_" + line.substring(1);
					seqSb = new StringBuilder();
				} else {
					seqSb.append(line);
				}
			}

			if (hostProName != null) {
				String seq = seqSb.toString();
				writer.println(hostProName);
				writer.println(seq);

				String cleavage = "KR";
				String shuffle = DecoyShuffler.generateShuffledDecoy(seq, cleavage, 42);
				writer.println(">shuf_" + hostProName.substring(1));
				writer.println(shuffle);

				writer.println(refHostProName);
				writer.println(new StringBuilder(seq).reverse().toString());

				writer.println(">rev_shuf_" + hostProName.substring(1));
				writer.println(new StringBuilder(shuffle).reverse().toString());
			}

			hostReader.close();

			writer.close();
		} catch (IOException e) {

		}
	}
	
	public static void main(String[] args) {
		FragpipePSMTsvReader.combineToFasta("Z:\\Test\\BellAshley\\MetaPilot\\fragpipe_MAG_search", 
				"C:\\Users\\daz25hak\\MetaProAI\\database\\uniprotkb_Human_AND_reviewed_true_AND_m_2025_09_01.fasta");
	}
}
