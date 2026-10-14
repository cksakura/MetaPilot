package uk.ac.quadram.metapilot.dbSearch.sage;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;

import org.apache.commons.io.FileUtils;

public class SagePsmReader {

	/**
	 * PSM-level FDR applied to the PSMs that reach MetaPilot's peptide/PSM layer -
	 * the same 1% Sage's own LFQ stage applies to peptide_q. Before this existed
	 * every PSM of a peptide listed in lfq.tsv was kept regardless of its own
	 * spectrum_q (on 12mix YS-194-1-1, 38,640 of 169,510 target PSMs with
	 * peptide_q <= 0.01 had spectrum_q > 0.01).
	 */
	public static final double PSM_Q_VALUE_THRESHOLD = 0.01;

	private SagePSM[] psms;

	private HashMap<String, Integer> ms2CountMap;

	public SagePsmReader(String in) {
		this(new File(in));
	}

	public SagePsmReader(File in) {
		read(in);
	}

	/**
	 * The spectra count in each raw file is unknown
	 * 
	 * @param in
	 */
	private void read(File in) {
		this.ms2CountMap = new HashMap<String, Integer>();
		ArrayList<SagePSM> list = new ArrayList<SagePSM>();
		BufferedReader reader = null;
		try {
			reader = new BufferedReader(new FileReader(in));
			String line = reader.readLine();

			String[] title = line.split("\t");

			int fileNameID = -1;
			int scanID = -1;
			int expMHID = -1;
			int chargeID = -1;
			int qValueID = -1;
			int sequenceID = -1;
			int calMHID = -1;
			int massPPMID = -1;
			int fragPPMID = -1;
			int hyperScoreID = -1;
			int disScoreID = -1;
			int specificitID = -1;
			int proteinID = -1;
			int missID = -1;
			int rtID = -1;
			int peptideQValueID = -1;
			int proteinQValueID = -1;
			int ms2IntensityID = -1;
			// Optional (QPX export / PSM filtering); absent columns leave defaults
			int labelID = -1;
			int rankID = -1;
			int posteriorErrorID = -1;
			int ionMobilityID = -1;

			for (int i = 0; i < title.length; i++) {
				if (title[i].equals("filename")) {
					fileNameID = i;
				} else if (title[i].equals("scannr")) {
					scanID = i;
				} else if (title[i].equals("expmass")) {
					expMHID = i;
				} else if (title[i].equals("charge")) {
					chargeID = i;
				} else if (title[i].equals("spectrum_q")) {
					qValueID = i;
				} else if (title[i].equals("peptide")) {
					sequenceID = i;
				} else if (title[i].equals("calcmass")) {
					calMHID = i;
				} else if (title[i].equals("precursor_ppm")) {
					massPPMID = i;
				} else if (title[i].equals("rt")) {
					rtID = i;
				} else if (title[i].equals("hyperscore")) {
					hyperScoreID = i;
				} else if (title[i].equals("sage_discriminant_score")) {
					disScoreID = i;
				} else if (title[i].equals("peptide_q")) {
					peptideQValueID = i;
				} else if (title[i].equals("semi_enzymatic")) {
					specificitID = i;
				} else if (title[i].equals("proteins")) {
					proteinID = i;
				} else if (title[i].equals("protein_q")) {
					proteinQValueID = i;
				} else if (title[i].equals("missed_cleavages")) {
					missID = i;
				} else if (title[i].equals("fragment_ppm")) {
					fragPPMID = i;
				} else if (title[i].equals("ms2_intensity")) {
					ms2IntensityID = i;
				} else if (title[i].equals("label")) {
					labelID = i;
				} else if (title[i].equals("rank")) {
					rankID = i;
				} else if (title[i].equals("posterior_error")) {
					posteriorErrorID = i;
				} else if (title[i].equals("ion_mobility")) {
					ionMobilityID = i;
				}
			}

			while ((line = reader.readLine()) != null) {
				String[] cs = line.split("\t");
				if (cs.length == title.length) {
					String fileName = cs[fileNameID];
					if (fileName.endsWith("mzML")) {
						fileName = fileName.substring(0, fileName.length() - ".mzML".length());
					}

					int scan = 0;
					int scanIndex = cs[scanID].indexOf("scan=");
					if (scanIndex > 0) {
						scan = Integer.parseInt(cs[scanID].substring(scanIndex + "scan=".length()));
					}

					String modseq = cs[sequenceID];
					StringBuilder sb = new StringBuilder();
					boolean mod = false;
					for (char aa : modseq.toCharArray()) {
						if (aa == '[') {
							mod = true;
						} else if (aa == ']') {
							mod = false;
						} else {
							if (!mod) {
								sb.append(aa);
							}
						}
					}
					double expMass = Double.parseDouble(cs[expMHID]);
					double calMass = Double.parseDouble(cs[calMHID]);
					double massPPM = Double.parseDouble(cs[massPPMID]);
					double fragPPM = Double.parseDouble(cs[fragPPMID]);
					int charge = Integer.parseInt(cs[chargeID]);
					int miss = Integer.parseInt(cs[missID]);
					int specificity = Integer.parseInt(cs[specificitID]);
					double rt = Double.parseDouble(cs[rtID]);
					double qValue = Double.parseDouble(cs[qValueID]);
					double hyperscore = Double.parseDouble(cs[hyperScoreID]);
					double disScore = Double.parseDouble(cs[disScoreID]);
					boolean isTarget = true;
					String[] proteins = cs[proteinID].split(";");
					// Sage's OWN decoys, which the identification output retains. Not
					// shuf_: those are targets as far as Sage is concerned, and they
					// carry genome-level FDR after quantification instead - see
					// SageDecoyTags.
					if (SageDecoyTags.anyEngineDecoy(proteins)) {
						isTarget = false;
					}
					// Sage's own label (1 target, -1 decoy) is authoritative when present.
					// On 12mix YS-194-1-1 it agreed with the accession rule on every row.
					if (labelID >= 0) {
						isTarget = "1".equals(cs[labelID].trim());
					}
					double peptideQ = Double.parseDouble(cs[peptideQValueID]);
					double proteinQ = Double.parseDouble(cs[proteinQValueID]);

					double ms2Intensity = Double.parseDouble(cs[ms2IntensityID]);

					SagePSM psm = new SagePSM(fileName, scan, sb.toString(), modseq, expMass, calMass, massPPM, fragPPM,
							charge, miss, specificity, rt, qValue, hyperscore, disScore, proteins, isTarget, peptideQ,
							proteinQ, ms2Intensity);
					if (rankID >= 0) {
						psm.setRank(parseIntOr(cs[rankID], 1));
					}
					if (posteriorErrorID >= 0) {
						psm.setPosteriorErrorLog10(parseDoubleOr(cs[posteriorErrorID], Double.NaN));
					}
					if (ionMobilityID >= 0) {
						psm.setIonMobility(parseDoubleOr(cs[ionMobilityID], Double.NaN));
					}
					list.add(psm);
				}
			}
			reader.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

		this.psms = list.toArray(new SagePSM[list.size()]);
	}

	public SagePSM[] getPsms() {
		return psms;
	}

	/**
	 * The PSMs MetaPilot treats as identifications: Sage targets, top-ranked
	 * candidate, spectrum_q <= {@link #PSM_Q_VALUE_THRESHOLD}.
	 */
	public static boolean isConfidentPsm(SagePSM psm) {
		return psm != null && psm.isTarget() && psm.getRank() <= 1
				&& psm.getqValue() <= PSM_Q_VALUE_THRESHOLD;
	}

	private static int parseIntOr(String s, int fallback) {
		try {
			return Integer.parseInt(s.trim());
		} catch (RuntimeException e) {
			return fallback;
		}
	}

	private static double parseDoubleOr(String s, double fallback) {
		try {
			return Double.parseDouble(s.trim());
		} catch (RuntimeException e) {
			return fallback;
		}
	}
	
	public SagePSM[] getPsms(double qvalueThreshold) {
		if (psms == null || psms.length == 0) {
			return new SagePSM[0];
		}
		ArrayList<SagePSM> filtered = new ArrayList<>();
		for (SagePSM psm : psms) {
			if (psm == null) {
				continue;
			}
			if (psm.getPeptideQ() < qvalueThreshold) {
				filtered.add(psm);
			}
		}
		return filtered.toArray(new SagePSM[filtered.size()]);
	}
	
	public SagePSM[] getPsms(double specQvalueThreshold, double pepQvalueThreshold) {
		if (psms == null || psms.length == 0) {
			return new SagePSM[0];
		}
		ArrayList<SagePSM> filtered = new ArrayList<>();
		for (SagePSM psm : psms) {
			if (psm == null) {
				continue;
			}
			if (psm.getqValue() < specQvalueThreshold && psm.getPeptideQ() < pepQvalueThreshold) {
				filtered.add(psm);
			}
		}
		return filtered.toArray(new SagePSM[filtered.size()]);
	}

	public HashMap<String, Integer> getMs2CountMap() {
		return ms2CountMap;
	}

	public static void filterFdr(String in, String out, double fdr) throws IOException {

		HashMap<String, float[]> scanNameMap = new HashMap<String, float[]>();
		BufferedReader reader = new BufferedReader(new FileReader(in));
		String line = reader.readLine();
		String[] title = line.split("\t");
		while ((line = reader.readLine()) != null) {

			String[] cs = line.split("\t");
			if (cs.length == title.length) {
				float qvalue = Float.parseFloat(cs[4]);
				if (cs[15].equals("target")) {
					scanNameMap.put(cs[0], new float[] { qvalue, 1.0f });
				} else {
					scanNameMap.put(cs[0], new float[] { qvalue, 0.0f });
				}
			}

		}
		reader.close();

		String[] scanNames = scanNameMap.keySet().toArray(new String[scanNameMap.size()]);
		Arrays.sort(scanNames, new Comparator<String>() {

			@Override
			public int compare(String o1, String o2) {
				// TODO Auto-generated method stub
				float[] f1 = scanNameMap.get(o1);
				float[] f2 = scanNameMap.get(o2);

				if (f1[0] < f2[0]) {
					return -1;
				} else if (f1[0] > f2[0]) {
					return 1;
				}

				return 0;
			}
		});

		int target = 0;
		int decoy = 0;
		boolean begin = false;
		for (int i = 0; i < scanNames.length; i++) {
			if (begin) {
				scanNameMap.remove(scanNames[i]);
			} else {
				float[] score = scanNameMap.get(scanNames[i]);
				if (score[1] == 0.0f) {
					decoy++;
				} else {
					target++;
				}

				if (target * fdr < decoy) {
					scanNameMap.remove(scanNames[i]);
					begin = true;
				}
			}
		}

		PrintWriter writer = new PrintWriter(out);
		reader = new BufferedReader(new FileReader(in));

		writer.println(reader.readLine());

		while ((line = reader.readLine()) != null) {
			String[] cs = line.split("\t");
			if (cs.length == title.length) {
				if (scanNameMap.containsKey(cs[0]) && cs[15].equals("target")) {
					writer.println(line);
				}
			}
		}
		reader.close();
		writer.close();
	}

	public static HashSet<String> getPepSet(String in) throws IOException {
		HashSet<String> set1 = new HashSet<String>();
		BufferedReader reader1 = new BufferedReader(new FileReader(in));
		String[] title = reader1.readLine().split("\t");
		String line = null;
		while ((line = reader1.readLine()) != null) {
			String[] cs = line.split("\t");
			if (cs.length == title.length) {
				set1.add(cs[5]);
			}
		}
		reader1.close();
		return set1;
	}
	
	public static void compare(String s1, String s2) {
		SagePsmReader reader1 = new SagePsmReader(s1);
		SagePSM[] psms1 = reader1.getPsms(0.01, 1.0);

		SagePsmReader reader2 = new SagePsmReader(s2);
		SagePSM[] psms2 = reader2.getPsms(0.01, 1.0);
		System.out.println(psms1.length + "\t" + psms2.length);

		HashMap<String, HashSet<String>> pepGenomeMap1 = new HashMap<String, HashSet<String>>();
		for (SagePSM psm : psms1) {
			if (psm.isTarget()) {
				String seq = psm.getBaseSequence();
				String[] proteins = psm.getProteins();
				for (String pro : proteins) {
					String genome = pro.substring(0, pro.indexOf("_"));
					pepGenomeMap1.computeIfAbsent(seq, k -> new HashSet<>()).add(genome);
				}
			}
		}

		HashMap<String, HashSet<String>> pepGenomeMap2 = new HashMap<String, HashSet<String>>();
		for (SagePSM psm : psms2) {
			if (psm.isTarget()) {
				String seq = psm.getBaseSequence();
				String[] proteins = psm.getProteins();
				for (String pro : proteins) {
					String genome = pro.substring(0, pro.indexOf("_"));
					pepGenomeMap2.computeIfAbsent(seq, k -> new HashSet<>()).add(genome);
				}
			}
		}
		HashSet<String> set = new HashSet<String>();
		set.addAll(pepGenomeMap1.keySet());
		set.addAll(pepGenomeMap2.keySet());
		int common = pepGenomeMap1.size() + pepGenomeMap2.size() - set.size();
		System.out.println(pepGenomeMap1.size() + "\t" + pepGenomeMap2.size());
		System.out.println((pepGenomeMap1.size() - common) + "\t" + (pepGenomeMap2.size() - common));

		HashMap<String, HashSet<String>> genomeMap = new HashMap<String, HashSet<String>>();
		for (String key : pepGenomeMap2.keySet()) {
			if (!pepGenomeMap1.containsKey(key)) {
				HashSet<String> genomeSet = pepGenomeMap2.get(key);
				for (String genome : genomeSet) {
					genomeMap.computeIfAbsent(genome, k -> new HashSet<>()).add(key);
				}
			}
		}
		for (String genome : genomeMap.keySet()) {
//			System.out.println(genome + "\t" + genomeMap.get(genome).size());
		}

		for (SagePSM psm : psms2) {
			String seq = psm.getBaseSequence();
			if (!pepGenomeMap1.containsKey(seq)) {
//				System.out.println(seq+"\t"+psm.getScan());
			}
		}
	}
	
	public static void fastaTest(String in) {
		try (BufferedReader reader = new BufferedReader(new FileReader(in))) {
			String line = null;
			HashMap<String, HashSet<String>> genomeProMap = new HashMap<String, HashSet<String>>();
			while ((line = reader.readLine()) != null) {
				if (line.startsWith(">")) {
					int id = line.indexOf(" ");
					if (id > 1) {
						String pro = line.substring(1, id);
						String genome = pro.substring(0, pro.indexOf("_"));
						genomeProMap.computeIfAbsent(genome, k -> new HashSet<>()).add(pro);
					} else {
						System.out.println(line);
					}
				}
			}
			reader.close();

			for (String genome : genomeProMap.keySet()) {
				System.out.println(genome + "\t" + genomeProMap.get(genome).size());
			}

		} catch (IOException e) {

		}
	}


	public static void main(String[] args) {
		// TODO Auto-generated method stub
		
//		SagePsmReader.fastaTest("Z:\\PRIDE_dataset\\PXD006118_3mix_type\\Run1\\MetaPilot\\sage_MAG_search_09\\old\\combined\\combined.fasta");
		
		File[] folders = (new File("Z:\\PRIDE_dataset\\PXD012724_rapidAim1\\MetaPilot\\sage_MAG_search_rib"))
				.listFiles();
		
		File folder = new File("Z:\\PRIDE_dataset\\PXD060335_campi2\\MetaPilot");
		File f1 = new File(folder, "sage_MAG_search");
		File f2 = new File(folder, "sage_MAG_search_rib");
		String fileName = "20230815_CAMPI2_CEB_F_M3_R1_QE3_2ug_1D_DirectLoad_RP180grad";
		
		File sub1 = new File(f1, fileName);
		File sub2 = new File(f2, fileName);
		
		File result1 = new File(sub1, "results.sage.tsv");
		File result2 = new File(sub2, "results.sage.tsv");
		
//		SagePsmReader.compare(result1.getAbsolutePath(), result2.getAbsolutePath());
		System.out.println();
		
		File hap1 = new File(sub1, "hap");
		File hap2 = new File(sub2, "hap");
		
		File hapResult1 = new File(hap1, "results.sage.tsv");
		File hapResult2 = new File(hap2, "results.sage.tsv");
//		SagePsmReader.compare(hapResult1.getAbsolutePath(), hapResult2.getAbsolutePath());
		
//		SagePsmReader.compare("Z:\\PRIDE_dataset\\PXD004039\\MetaPilot\\sage_MAG_search\\Human_3\\hap\\marker4_highpep\\results.sage.tsv",
//				"Z:\\PRIDE_dataset\\PXD004039\\MetaPilot\\sage_MAG_search\\Human_3\\hap\\marker4\\results.sage.tsv"
//				);	
		System.out.println();
//		SagePsmReader.compare(
//				"Z:\\PRIDE_dataset\\PXD004039\\MetaPilot\\sage_MAG_search\\Human_3\\hap\\rib\\results.sage.tsv",
//				"Z:\\PRIDE_dataset\\PXD004039\\MetaPilot\\sage_MAG_search\\Human_3\\hap\\marker4\\results.sage.tsv");	
		
		for (int i = 0; i < folders.length; i++) {
			File resultFolder = new File(folders[i], "results.sage.tsv");
			if (resultFolder.exists()) {
				SagePsmReader reader = new SagePsmReader(resultFolder);
				System.out.println(folders[i].getName() + "\t" + reader.getPsms(0.01).length);
//				FileUtils.deleteQuietly(resultFolder);
			}
			
			File hapFile = new File(new File(folders[i], "hap"), "results.sage.tsv");
			if (hapFile.exists()) {
				SagePsmReader reader = new SagePsmReader(hapFile);
				System.out.println(folders[i].getName() + "\t" + reader.getPsms(0.01).length);
//				FileUtils.deleteQuietly(resultFolder);
			}
		}		
		
		
		for (int i = 0; i < folders.length; i++) {
			File resultFolder = new File(folders[i], "hap");
			if (resultFolder.exists()) {
				File fastaFile = new File(resultFolder, folders[i].getName() + ".fasta");
				if (fastaFile.exists()) {
//					FileUtils.deleteQuietly(fastaFile);
					/*
					try (BufferedReader reader = new BufferedReader(new FileReader(fastaFile))) {
						String line = null;
						int count = 0;
						HashSet<String> set = new HashSet<String>();
						while ((line = reader.readLine()) != null) {
							if (line.startsWith(">")) {
								count++;
								String genome = line.substring(1, line.indexOf("_"));
								set.add(genome);
							}
						}
						System.out.println(folders[i].getName() + "\t" + count + "\t" + set.size());
					} catch (IOException e) {

					}
					*/
				}
			}
		}
	}

}
