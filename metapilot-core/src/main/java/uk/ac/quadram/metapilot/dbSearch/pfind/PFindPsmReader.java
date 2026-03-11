/**
 * 
 */
package uk.ac.quadram.metapilot.dbSearch.pfind;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileFilter;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;

import org.apache.commons.io.FileUtils;

import uk.ac.quadram.metapilot.core.mod.PostTransModification;

/**
 * Use this for the reading of pFind.spectra, use PFindResultReader for the reading of pFind-filtered.spectra
 * 
 * @see bmi.med.uOttawa.metalab.dbSearch.pfind.PFindResultReader
 * @author Kai Cheng
 */
public class PFindPsmReader {

	private PFindPSM[] psms;

	private static final String REV = "REV_";

	private HashMap<String, Integer> ms2CountMap;
	private HashMap<String, PostTransModification> modMap;

	public PFindPsmReader(String in) {
		this(new File(in));
	}

	public PFindPsmReader(File in) {
		read(in, new HashMap<String, PostTransModification>());
	}

	public PFindPsmReader(String in, HashMap<String, PostTransModification> modMap) {
		this(new File(in), modMap);
	}

	public PFindPsmReader(File in, HashMap<String, PostTransModification> modMap) {
		read(in, modMap);
	}

	public PFindPsmReader(String in, String ms2SpecDir) {
		this(new File(in), new File(ms2SpecDir));
	}

	public PFindPsmReader(File in, File ms2SpecDir) {
		read(in, ms2SpecDir, new HashMap<String, PostTransModification>());
	}

	public PFindPsmReader(File in, File ms2SpecDir, HashMap<String, PostTransModification> modMap) {
		read(in, ms2SpecDir, modMap);
	}

	public PFindPsmReader(String in, String ms2SpecDir, HashMap<String, HashSet<String>> promap) {
		this(new File(in), new File(ms2SpecDir), promap, new HashMap<String, PostTransModification>());
	}

	public PFindPsmReader(File in, File ms2SpecDir, HashMap<String, HashSet<String>> promap,
			HashMap<String, PostTransModification> modMap) {
		try {
			read(in, ms2SpecDir, promap, modMap);
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

	/**
	 * The spectra count in each raw file is unknown
	 * 
	 * @param in
	 */
	private void read(File in, HashMap<String, PostTransModification> modMap) {
		this.ms2CountMap = new HashMap<String, Integer>();
		ArrayList<PFindPSM> list = new ArrayList<PFindPSM>();
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
			int massShiftID = -1;
			int rawScoreID = -1;
			int finalScoreID = -1;
			int modID = -1;
			int specificitID = -1;
			int proteinID = -1;
			int isTargetID = -1;
			int missID = -1;
			int fragMassShiftID = -1;

			for (int i = 0; i < title.length; i++) {
				if (title[i].equals("File_Name")) {
					fileNameID = i;
				} else if (title[i].equals("Scan_No")) {
					scanID = i;
				} else if (title[i].equals("Exp.MH+")) {
					expMHID = i;
				} else if (title[i].equals("Charge")) {
					chargeID = i;
				} else if (title[i].equals("Q-value")) {
					qValueID = i;
				} else if (title[i].equals("Sequence")) {
					sequenceID = i;
				} else if (title[i].equals("Calc.MH+")) {
					calMHID = i;
				} else if (title[i].equals("Mass_Shift(Exp.-Calc.)")) {
					massShiftID = i;
				} else if (title[i].equals("Raw_Score")) {
					rawScoreID = i;
				} else if (title[i].equals("Final_Score")) {
					finalScoreID = i;
				} else if (title[i].equals("Modification")) {
					modID = i;
				} else if (title[i].equals("Specificity")) {
					specificitID = i;
				} else if (title[i].equals("Proteins")) {
					proteinID = i;
				} else if (title[i].equals("Target/Decoy")) {
					isTargetID = i;
				} else if (title[i].equals("Miss.Clv.Sites")) {
					missID = i;
				} else if (title[i].equals("Avg.Frag.Mass.Shift")) {
					fragMassShiftID = i;
				}
			}

			while ((line = reader.readLine()) != null) {
				String[] cs = line.split("\t");
				if (cs.length == title.length) {

					String fileName = cs[fileNameID];
					if (fileName.indexOf(".") >= 0) {
						fileName = fileName.substring(0, fileName.indexOf("."));
					}

					int scan = Integer.parseInt(cs[scanID]);
					double expMH = Double.parseDouble(cs[expMHID]);
					int charge = Integer.parseInt(cs[chargeID]);
					double qValue = Double.parseDouble(cs[qValueID]);
					String sequence = cs[sequenceID];
					if (sequence.length() < 7) {
						continue;
					}
					if (qValue > 0.01) {
						continue;
					}

					double calMH = Double.parseDouble(cs[calMHID]);
					double massShift = Double.parseDouble(cs[massShiftID]);
					double rawScore = Double.parseDouble(cs[rawScoreID]);
					double finalScore = -Math.log10(Double.parseDouble(cs[finalScoreID]));
					String mod = cs[modID];
					int specificity = Integer.parseInt(cs[specificitID]);
					boolean isTarget = cs[isTargetID].equals("target");
					String[] proteins = cs[proteinID].split("/");

					// for short peptides (generally 6 aminoacids), they can belong to both target and decoy proteins
					ArrayList<String> decoyList = new ArrayList<String>();
					ArrayList<String> targetList = new ArrayList<String>();
					for (int i = 0; i < proteins.length; i++) {
						if (proteins[i].startsWith(REV)) {
							decoyList.add(proteins[i]);
						} else {
							targetList.add(proteins[i]);
						}
					}
					if (decoyList.size() > 0) {
						isTarget = false;
						proteins = decoyList.toArray(new String[decoyList.size()]);
					} else if (isTarget) {
						if (targetList.size() == 0) {
							continue;
						}
						proteins = targetList.toArray(new String[targetList.size()]);
					}

					int miss = Integer.parseInt(cs[missID]);
					double fragMassShift = Double.parseDouble(cs[fragMassShiftID]);

					PFindPSM psm = new PFindPSM(fileName, scan, expMH, charge, qValue, sequence, calMH, massShift,
							rawScore, finalScore, mod, modMap.get(mod), specificity, proteins, isTarget, miss,
							fragMassShift);
					list.add(psm);
				}
			}
			reader.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

		this.psms = list.toArray(new PFindPSM[list.size()]);
	}

	/**
	 * Get the information of retention time from the mgf files and add it to the
	 * psm
	 * 
	 * @param in
	 * @param rawDir
	 */
	private void read(File in, File rawDir, HashMap<String, PostTransModification> modMap) {

		this.ms2CountMap = new HashMap<String, Integer>();
		HashMap<String, HashMap<Integer, Double>> rtmap = new HashMap<String, HashMap<Integer, Double>>();

		File rtFile = new File(rawDir, "rt.txt");
		if (rtFile.exists() && rtFile.length() > 0) {
			BufferedReader reader = null;
			try {
				reader = new BufferedReader(new FileReader(rtFile));
				String line = reader.readLine();
				String[] title = line.split("\t");
				while ((line = reader.readLine()) != null) {
					String[] cs = line.split("\t");
					if (cs.length == title.length) {
						if (rtmap.containsKey(cs[0])) {
							rtmap.get(cs[0]).put(Integer.parseInt(cs[1]), Double.parseDouble(cs[2]));
						} else {
							HashMap<Integer, Double> map = new HashMap<Integer, Double>();
							map.put(Integer.parseInt(cs[1]), Double.parseDouble(cs[2]));
							rtmap.put(cs[0], map);
						}
					}
				}
				reader.close();

			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}

		} else {
			File[] ms2Files = rawDir.listFiles(new FileFilter() {

				@Override
				public boolean accept(File arg0) {
					// TODO Auto-generated method stub
					if (arg0.getName().endsWith("ms2")) {
						return true;
					}
					return false;
				}

			});

			for (int i = 0; i < ms2Files.length; i++) {

				String line = null;
				int scan = -1;
				double rt = -1;
				try {
					int ms2Count = 0;
					HashMap<Integer, Double> map = new HashMap<Integer, Double>();
					BufferedReader reader = new BufferedReader(new FileReader(ms2Files[i]));
					while ((line = reader.readLine()) != null) {
						String[] cs = line.split("\t");
						if (cs.length == 3) {

							if (cs[1].equals("RetTime")) {
								rt = Double.parseDouble(cs[2]);
								map.put(scan, rt / 60.0);
							}
						} else if (cs.length == 4) {
							if (cs[0].equals("S")) {
								scan = Integer.parseInt(cs[1]);
								ms2Count++;
							}
						}
					}
					reader.close();

					String name = ms2Files[i].getName();
					name = name.substring(0, name.lastIndexOf("."));

					rtmap.put(name, map);
					this.ms2CountMap.put(name, ms2Count);

				} catch (IOException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
			}
		}

		ArrayList<PFindPSM> list = new ArrayList<PFindPSM>();
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
			int massShiftID = -1;
			int rawScoreID = -1;
			int finalScoreID = -1;
			int modID = -1;
			int specificitID = -1;
			int proteinID = -1;
			int isTargetID = -1;
			int missID = -1;
			int fragMassShiftID = -1;

			for (int i = 0; i < title.length; i++) {
				if (title[i].equals("File_Name")) {
					fileNameID = i;
				} else if (title[i].equals("Scan_No")) {
					scanID = i;
				} else if (title[i].equals("Exp.MH+")) {
					expMHID = i;
				} else if (title[i].equals("Charge")) {
					chargeID = i;
				} else if (title[i].equals("Q-value")) {
					qValueID = i;
				} else if (title[i].equals("Sequence")) {
					sequenceID = i;
				} else if (title[i].equals("Calc.MH+")) {
					calMHID = i;
				} else if (title[i].equals("Mass_Shift(Exp.-Calc.)")) {
					massShiftID = i;
				} else if (title[i].equals("Raw_Score")) {
					rawScoreID = i;
				} else if (title[i].equals("Final_Score")) {
					finalScoreID = i;
				} else if (title[i].equals("Modification")) {
					modID = i;
				} else if (title[i].equals("Specificity")) {
					specificitID = i;
				} else if (title[i].equals("Proteins")) {
					proteinID = i;
				} else if (title[i].equals("Target/Decoy")) {
					isTargetID = i;
				} else if (title[i].equals("Miss.Clv.Sites")) {
					missID = i;
				} else if (title[i].equals("Avg.Frag.Mass.Shift")) {
					fragMassShiftID = i;
				}
			}

			while ((line = reader.readLine()) != null) {
				String[] cs = line.split("\t");
				if (cs.length == title.length) {

					String fileName = cs[fileNameID];
					fileName = fileName.substring(0, fileName.indexOf("."));
					int scan = Integer.parseInt(cs[scanID]);
					double expMH = Double.parseDouble(cs[expMHID]);
					int charge = Integer.parseInt(cs[chargeID]);
					double qValue = Double.parseDouble(cs[qValueID]);
					String sequence = cs[sequenceID];
					double calMH = Double.parseDouble(cs[calMHID]);
					double massShift = Double.parseDouble(cs[massShiftID]);
					double rawScore = Double.parseDouble(cs[rawScoreID]);
					double finalScore = Double.parseDouble(cs[finalScoreID]);
					String mod = cs[modID];
					int specificity = Integer.parseInt(cs[specificitID]);
					boolean isTarget = cs[isTargetID].equals("target");
					String[] proteins = cs[proteinID].split("/");

					// for short peptides (generally 6 aminoacids), they can belong to both target
					// and decoy proteins
					if (isTarget) {
						ArrayList<String> proList = new ArrayList<String>();
						for (int i = 0; i < proteins.length; i++) {
							if (!proteins[i].startsWith(REV)) {
								proList.add(proteins[i]);
							}
						}
						if (proList.size() == 0) {
							continue;
						}
						proteins = proList.toArray(new String[proList.size()]);
					}

					int miss = Integer.parseInt(cs[missID]);
					double fragMassShift = Double.parseDouble(cs[fragMassShiftID]);

					PFindPSM psm = new PFindPSM(fileName, scan, expMH, charge, qValue, sequence, calMH, massShift,
							rawScore, finalScore, mod, modMap.get(mod), specificity, proteins, isTarget, miss,
							fragMassShift);

					if (rtmap.containsKey(fileName)) {
						HashMap<Integer, Double> map = rtmap.get(fileName);
						if (map.containsKey(scan)) {
							psm.setRt(map.get(scan));
						}
					}

					list.add(psm);
				}
			}
			reader.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

		this.psms = list.toArray(new PFindPSM[list.size()]);
	}

	/**
	 * The protein information is not read from the "pFind-Filtered.spectra" file
	 * directly, but from the "pFind.protein" file.
	 * 
	 * @param in
	 * @param rawDir
	 * @param promap
	 * @throws IOException 
	 */
	private void read(File in, File ms2SpecDir, HashMap<String, HashSet<String>> promap,
			HashMap<String, PostTransModification> modMap) throws IOException {

		File ms2CountFile = new File(ms2SpecDir, "ms2Count.txt");
		this.ms2CountMap = new HashMap<String, Integer>();
		if (ms2CountFile.exists()) {
			BufferedReader reader = new BufferedReader(new FileReader(ms2CountFile));
			String line = reader.readLine();
			while ((line = reader.readLine()) != null) {
				String[] cs = line.split("\t");
				ms2CountMap.put(cs[0], Integer.parseInt(cs[1]));
			}
			reader.close();
		}

		File rtFile = new File(ms2SpecDir, "rt.txt");
		HashMap<String, HashMap<Integer, Double>> rtmap = new HashMap<String, HashMap<Integer, Double>>();
		if (rtFile.exists()) {
			BufferedReader reader = new BufferedReader(new FileReader(rtFile));
			String line = reader.readLine();
			while ((line = reader.readLine()) != null) {
				String[] cs = line.split("\t");
				HashMap<Integer, Double> rtmapi;
				if (rtmap.containsKey(cs[0])) {
					rtmapi = rtmap.get(cs[0]);
				} else {
					rtmapi = new HashMap<Integer, Double>();
					rtmap.put(cs[0], rtmapi);
				}
				rtmapi.put(Integer.parseInt(cs[1]), Double.parseDouble(cs[2]));
			}
			reader.close();

		} else {

			File[] ms2Files = ms2SpecDir.listFiles(new FileFilter() {

				@Override
				public boolean accept(File arg0) {
					// TODO Auto-generated method stub
					if (arg0.getName().endsWith("ms2")) {
						return true;
					}
					return false;
				}

			});

			PrintWriter rtWriter = new PrintWriter(rtFile);
			PrintWriter ms2CountWriter = new PrintWriter(ms2CountFile);

			rtWriter.println("FileName\tScannum\tRt");
			ms2CountWriter.println("FileName\tMS2 count");

			for (int i = 0; i < ms2Files.length; i++) {

				String line = null;
				int scan = -1;
				double rt = -1;
				int ms2Count = 0;
				HashMap<Integer, Double> map = new HashMap<Integer, Double>();
				BufferedReader reader = new BufferedReader(new FileReader(ms2Files[i]));
				while ((line = reader.readLine()) != null) {
					String[] cs = line.split("\t");
					if (cs.length == 3) {

						if (cs[1].equals("RetTime")) {
							rt = Double.parseDouble(cs[2]);
							map.put(scan, rt / 60.0);
						}
					} else if (cs.length == 4) {
						if (cs[0].equals("S")) {
							scan = Integer.parseInt(cs[1]);
							ms2Count++;
						}
					}
				}
				reader.close();

				String name = ms2Files[i].getName();
				name = name.substring(0, name.lastIndexOf("."));

				rtmap.put(name, map);
				this.ms2CountMap.put(name, ms2Count);

				for (Integer scannumi : map.keySet()) {
					rtWriter.println(name + "\t" + scannumi + "\t" + map.get(scannumi));
				}
				ms2CountWriter.println(name + "\t" + ms2Count);
			}

			rtWriter.close();
			ms2CountWriter.close();
		}

		ArrayList<PFindPSM> list = new ArrayList<PFindPSM>();
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
			int massShiftID = -1;
			int rawScoreID = -1;
			int finalScoreID = -1;
			int modID = -1;
			int specificitID = -1;
//			int proteinID = -1;
			int isTargetID = -1;
			int missID = -1;
			int fragMassShiftID = -1;

			for (int i = 0; i < title.length; i++) {
				if (title[i].equals("File_Name")) {
					fileNameID = i;
				} else if (title[i].equals("Scan_No")) {
					scanID = i;
				} else if (title[i].equals("Exp.MH+")) {
					expMHID = i;
				} else if (title[i].equals("Charge")) {
					chargeID = i;
				} else if (title[i].equals("Q-value")) {
					qValueID = i;
				} else if (title[i].equals("Sequence")) {
					sequenceID = i;
				} else if (title[i].equals("Calc.MH+")) {
					calMHID = i;
				} else if (title[i].equals("Mass_Shift(Exp.-Calc.)")) {
					massShiftID = i;
				} else if (title[i].equals("Raw_Score")) {
					rawScoreID = i;
				} else if (title[i].equals("Final_Score")) {
					finalScoreID = i;
				} else if (title[i].equals("Modification")) {
					modID = i;
				} else if (title[i].equals("Specificity")) {
					specificitID = i;
				} else if (title[i].equals("Proteins")) {
//					proteinID = i;
				} else if (title[i].equals("Target/Decoy")) {
					isTargetID = i;
				} else if (title[i].equals("Miss.Clv.Sites")) {
					missID = i;
				} else if (title[i].equals("Avg.Frag.Mass.Shift")) {
					fragMassShiftID = i;
				}
			}

			while ((line = reader.readLine()) != null) {
				String[] cs = line.split("\t");
				if (cs.length == title.length) {

					String fileName = cs[fileNameID];
					fileName = fileName.substring(0, fileName.indexOf("."));
					int scan = Integer.parseInt(cs[scanID]);
					double expMH = Double.parseDouble(cs[expMHID]);
					int charge = Integer.parseInt(cs[chargeID]);
					double qValue = Double.parseDouble(cs[qValueID]);
					String sequence = cs[sequenceID];
					double calMH = Double.parseDouble(cs[calMHID]);
					double massShift = Double.parseDouble(cs[massShiftID]);
					double rawScore = Double.parseDouble(cs[rawScoreID]);
					double finalScore = Double.parseDouble(cs[finalScoreID]);
					String mod = cs[modID];
					int specificity = Integer.parseInt(cs[specificitID]);

					String[] proteins;

					if (promap.containsKey(sequence)) {
						HashSet<String> plist = promap.get(sequence);
						proteins = plist.toArray(new String[plist.size()]);
					} else {
						/*
						 * String[] proteins = cs[proteinID].split("/"); for (int i = 0; i <
						 * proteins.length; i++) { prosb.append(proteins[i]).append(";"); }
						 */

						continue;
					}

					boolean isTarget = cs[isTargetID].equals("target");
					int miss = Integer.parseInt(cs[missID]);
					double fragMassShift = Double.parseDouble(cs[fragMassShiftID]);

					PFindPSM psm = new PFindPSM(fileName, scan, expMH, charge, qValue, sequence, calMH, massShift,
							rawScore, finalScore, mod, modMap.get(mod), specificity, proteins, isTarget, miss,
							fragMassShift);

					if (rtmap.containsKey(fileName)) {
						HashMap<Integer, Double> map = rtmap.get(fileName);
						if (map.containsKey(scan)) {
							psm.setRt(map.get(scan));
						}
					}

					list.add(psm);
				}
			}
			reader.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

		this.psms = list.toArray(new PFindPSM[list.size()]);
	}

	/**
	 * The protein information is not read from the "pFind-Filtered.spectra" file
	 * directly, but from the "pFind.protein" file.
	 * 
	 * @param in
	 * @param rawDir
	 * @param promap
	 */
	@SuppressWarnings("unused")
	private void read(File in, HashMap<String, HashSet<String>> promap,
			HashMap<String, Integer> ms2CountMap, HashMap<String, HashMap<Integer, Double>> rtmap,
			HashMap<String, PostTransModification> modMap) {

		this.ms2CountMap = ms2CountMap;
		this.modMap = modMap;

		ArrayList<PFindPSM> list = new ArrayList<PFindPSM>();
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
			int massShiftID = -1;
			int rawScoreID = -1;
			int finalScoreID = -1;
			int modID = -1;
			int specificitID = -1;
			int proteinID = -1;
			int isTargetID = -1;
			int missID = -1;
			int fragMassShiftID = -1;

			for (int i = 0; i < title.length; i++) {
				if (title[i].equals("File_Name")) {
					fileNameID = i;
				} else if (title[i].equals("Scan_No")) {
					scanID = i;
				} else if (title[i].equals("Exp.MH+")) {
					expMHID = i;
				} else if (title[i].equals("Charge")) {
					chargeID = i;
				} else if (title[i].equals("Q-value")) {
					qValueID = i;
				} else if (title[i].equals("Sequence")) {
					sequenceID = i;
				} else if (title[i].equals("Calc.MH+")) {
					calMHID = i;
				} else if (title[i].equals("Mass_Shift(Exp.-Calc.)")) {
					massShiftID = i;
				} else if (title[i].equals("Raw_Score")) {
					rawScoreID = i;
				} else if (title[i].equals("Final_Score")) {
					finalScoreID = i;
				} else if (title[i].equals("Modification")) {
					modID = i;
				} else if (title[i].equals("Specificity")) {
					specificitID = i;
				} else if (title[i].equals("Proteins")) {
					proteinID = i;
				} else if (title[i].equals("Target/Decoy")) {
					isTargetID = i;
				} else if (title[i].equals("Miss.Clv.Sites")) {
					missID = i;
				} else if (title[i].equals("Avg.Frag.Mass.Shift")) {
					fragMassShiftID = i;
				}
			}

			while ((line = reader.readLine()) != null) {
				String[] cs = line.split("\t");
				if (cs.length == title.length) {

					String fileName = cs[fileNameID];
					fileName = fileName.substring(0, fileName.indexOf("."));
					int scan = Integer.parseInt(cs[scanID]);
					double expMH = Double.parseDouble(cs[expMHID]);
					int charge = Integer.parseInt(cs[chargeID]);
					double qValue = Double.parseDouble(cs[qValueID]);
					String sequence = cs[sequenceID];
					double calMH = Double.parseDouble(cs[calMHID]);
					double massShift = Double.parseDouble(cs[massShiftID]);
					double rawScore = Double.parseDouble(cs[rawScoreID]);
					double finalScore = Double.parseDouble(cs[finalScoreID]);
					String mod = cs[modID];
					int specificity = Integer.parseInt(cs[specificitID]);

					String[] proteins;

					if (promap.containsKey(sequence)) {
						HashSet<String> plist = promap.get(sequence);
						proteins = plist.toArray(new String[plist.size()]);
					} else {
						/*
						 * String[] proteins = cs[proteinID].split("/"); for (int i = 0; i <
						 * proteins.length; i++) { prosb.append(proteins[i]).append(";"); }
						 */

						continue;
					}

					boolean isTarget = cs[isTargetID].equals("target");
					int miss = Integer.parseInt(cs[missID]);
					double fragMassShift = Double.parseDouble(cs[fragMassShiftID]);

					PFindPSM psm = new PFindPSM(fileName, scan, expMH, charge, qValue, sequence, calMH, massShift,
							rawScore, finalScore, mod, modMap.get(mod), specificity, proteins, isTarget, miss,
							fragMassShift);

					if (rtmap.containsKey(fileName)) {
						HashMap<Integer, Double> map = rtmap.get(fileName);
						if (map.containsKey(scan)) {
							psm.setRt(map.get(scan));
						}
					}

					list.add(psm);
				}
			}
			reader.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

		this.psms = list.toArray(new PFindPSM[list.size()]);
	}

	public PFindPSM[] getPsms() {
		return psms;
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
	
	public HashMap<String, PostTransModification> getModMap() {
		return modMap;
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
	
	public static void compareFolder(String s1, String s2) {
		HashMap<String, String> map1 = new HashMap<String, String>();
		File[] files1 = (new File(s1)).listFiles();
		for (int i = 0; i < files1.length; i++) {
			File resultFolder = new File(files1[i], "hap");
			File resultFile = new File(resultFolder, "pFind-Filtered.spectra");
			if (resultFile.exists()) {
				map1.put(files1[i].getName(), resultFile.getAbsolutePath());
			}
		}
		System.out.println(map1.size());

		HashMap<String, String> map2 = new HashMap<String, String>();
		File[] files2 = (new File(s2)).listFiles();
		for (int i = 0; i < files2.length; i++) {
			File resultFolder = new File(files2[i], "hap");
			File resultFile = new File(resultFolder, "pFind-Filtered.spectra");
			if (resultFile.exists()) {
				map2.put(files2[i].getName(), resultFile.getAbsolutePath());
			}
		}
		System.out.println(map2.size());

		for (String key : map1.keySet()) {
			if (map2.containsKey(key)) {
				System.out.println(key);
				compare(map1.get(key), map2.get(key));
				System.out.println();

				File fasta1 = new File(s1 + "\\" + key + "\\hap\\" + key + ".fasta");
				File fasta2 = new File(s2 + "\\" + key + "\\hap\\" + key + ".fasta");

				compareFasta(fasta1.getAbsolutePath(), fasta2.getAbsolutePath());
				System.out.println();
				
				File result1 = new File(s1 + "\\" + key + "\\pFind-Filtered.spectra");
				File result2 = new File(s2 + "\\" + key + "\\pFind-Filtered.spectra");
				
				compare(result1.getAbsolutePath(), result2.getAbsolutePath());
				System.out.println();
			}
		}
		System.out.println();
	}
	
	public static void compareFasta(String s1, String s2) {
		try {
			HashSet<String> genomeSet1 = new HashSet<String>();
			HashSet<String> proSet1 = new HashSet<String>();
			BufferedReader reader1 = new BufferedReader(new FileReader(s1));
			String line = null;
			while ((line = reader1.readLine()) != null) {
				if (line.startsWith(">MGYG")) {
					String pro = line.substring(1, line.indexOf(" "));
					String genome = pro.substring(0, pro.indexOf("_"));
					genomeSet1.add(genome);
					proSet1.add(pro);
				}
			}
			reader1.close();

			HashSet<String> genomeSet2 = new HashSet<String>();
			HashSet<String> proSet2 = new HashSet<String>();
			BufferedReader reader2 = new BufferedReader(new FileReader(s2));
			while ((line = reader2.readLine()) != null) {
				if (line.startsWith(">MGYG")) {
					String pro = line.substring(1, line.indexOf(" "));
					String genome = pro.substring(0, pro.indexOf("_"));
					genomeSet2.add(genome);
					proSet2.add(pro);
				}
			}
			reader2.close();

			HashSet<String> allGenomeSet = new HashSet<String>();
			HashSet<String> allProSet = new HashSet<String>();
			allGenomeSet.addAll(genomeSet1);
			allGenomeSet.addAll(genomeSet2);
			allProSet.addAll(proSet1);
			allProSet.addAll(proSet2);

			System.out.println(genomeSet1.size() + "\t" + genomeSet2.size() + "\t" + allGenomeSet.size());
			System.out.println(proSet1.size() + "\t" + proSet2.size() + "\t" + allProSet.size());

		} catch (IOException e) {

		}
	}

	public static void compare(String s1, String s2) {
		PFindPsmReader reader1 = new PFindPsmReader(s1);
		PFindPSM[] psms1 = reader1.getPsms();

		PFindPsmReader reader2 = new PFindPsmReader(s2);
		PFindPSM[] psms2 = reader2.getPsms();
		System.out.println(psms1.length + "\t" + psms2.length);

		HashMap<String, HashSet<String>> pepGenomeMap1 = new HashMap<String, HashSet<String>>();
		for (PFindPSM psm : psms1) {
			if (psm.isTarget()) {
				String seq = psm.getSequence();
				String[] proteins = psm.getProteins();
				for (String pro : proteins) {
					String genome = pro.substring(0, pro.indexOf("_"));
					pepGenomeMap1.computeIfAbsent(seq, k -> new HashSet<>()).add(genome);
				}
			}
		}

		HashMap<String, HashSet<String>> pepGenomeMap2 = new HashMap<String, HashSet<String>>();
		for (PFindPSM psm : psms2) {
			if (psm.isTarget()) {
				String seq = psm.getSequence();
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

		for (PFindPSM psm : psms2) {
			String seq = psm.getSequence();
			if (!pepGenomeMap1.containsKey(seq)) {
//				System.out.println(seq+"\t"+psm.getScan());
			}
		}
	}
	
	public static void cleanFolder(String in) {
		File[] files = (new File(in)).listFiles();
		for (int i = 0; i < files.length; i++) {
			File hapFile = new File(files[i], "hap");
			if (hapFile.exists()) {
				File[] subFiles = hapFile.listFiles();
				for (int j = 0; j < subFiles.length; j++) {
					if (!subFiles[j].getName().startsWith("marker_genome")) {
						FileUtils.deleteQuietly(subFiles[j]);
					}
					
					if(subFiles[j].getName().endsWith("fasta")) {
						FileUtils.deleteQuietly(subFiles[j]);
					}
				}

//				FileUtils.deleteQuietly(new File(files[i], "pFind-Filtered.spectra"));
//				FileUtils.deleteQuietly(new File(files[i], "pFind.spectra"));
//				File fastaFile = new File(hapFile, files[i].getName() + ".fasta");
//				FileUtils.deleteQuietly(fastaFile);
			} else {
				FileUtils.deleteQuietly(files[i]);
			}
		}
	}
	
	public static void copyCfg(String in) {
		File[] files = (new File(in)).listFiles();
		for (int i = 0; i < files.length; i++) {
			File hapFile = new File(files[i], "hap");
			if (hapFile.exists()) {
				File marker4HapFile = new File(files[i], "marker4_85");
//				marker4HapFile.mkdir();

				try {
					BufferedReader reader = new BufferedReader(new FileReader(new File(hapFile, "pFind.cfg")));
					PrintWriter writer = new PrintWriter(new File(marker4HapFile, "pFind.cfg"));
//					File fastaFile = new File("C:\\Users\\daz25hak\\MetaProAI\\database\\MAG\\human-gut\\v2.0.2\\4markers_highpep.fasta");
					File fastaFile = new File(marker4HapFile, files[i].getName()+".fasta");
					String line = null;
					while ((line = reader.readLine()) != null) {
						if (line.startsWith("fastapath")) {
							writer.println("fastapath=" + fastaFile.getAbsolutePath());
						} else if (line.startsWith("outputpath")) {
							writer.println("outputpath=" + files[i].getAbsolutePath());
						} else if (line.startsWith("open=")) {
							writer.println("open=1");
						} else {
							writer.println(line);
						}
					}
					reader.close();
					writer.close();
					
					writer = new PrintWriter(new File(marker4HapFile, "pFind.bat"));
					writer.println("C:");
					writer.println("cd C:\\Users\\daz25hak\\MetaProAI\\Resources\\DBReducer\\pFind3\\bin");
					writer.println("searcher.exe "+(new File(marker4HapFile, "pFind.cfg")).getAbsolutePath());
					writer.close();
					
				} catch (IOException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
			}
		}
	}
	
	public static void copyCfgFasta(String in, String fastaDb) {
		File[] files = (new File(in)).listFiles();
		for (int i = 0; i < files.length; i++) {
			File marker4HapFile = new File(files[i], "hap");
			if (marker4HapFile.exists()) {
				
				try {

					BufferedReader selectReader = new BufferedReader(
							new FileReader(new File(marker4HapFile, "strict_selected.tsv")));
					File fastaFile = new File(marker4HapFile, files[i].getName() + ".fasta");
//					if (!fastaFile.exists()) {
						PrintWriter selectWriter = new PrintWriter(fastaFile);

						HashSet<String> genomeSet = new HashSet<>();
						String line = selectReader.readLine(); // header
						while ((line = selectReader.readLine()) != null) {
							line = line.trim();
							if (line.isEmpty() || line.startsWith("#")) {
								break;
							}
							String[] parts = line.split("\t");
							if (parts.length > 0 && !parts[0].trim().isEmpty()) {
								genomeSet.add(parts[0].trim());
							}
						}
						selectReader.close();
						System.out.println(files[i].getName() + "\t" + genomeSet.size());

						for (String genome : genomeSet) {
							File genomeDbFile = new File(fastaDb, genome + ".faa");
							if (!genomeDbFile.exists()) {
								genomeDbFile = new File(fastaDb, genome + ".fasta");
							}
							if (genomeDbFile.exists()) {
								BufferedReader reader = new BufferedReader(new FileReader(genomeDbFile));
								while ((line = reader.readLine()) != null) {
									selectWriter.println(line);
								}
								reader.close();
							}
						}
						selectWriter.close();

						BufferedReader reader = new BufferedReader(
								new FileReader(new File(marker4HapFile, "pFind.cfg")));
						PrintWriter writer = new PrintWriter(new File(files[i], "pFind.cfg"));
						while ((line = reader.readLine()) != null) {
							if (line.startsWith("fastapath")) {
								writer.println("fastapath=" + fastaFile.getAbsolutePath());
							} else if (line.startsWith("outputpath")) {
								writer.println("outputpath=" + files[i].getAbsolutePath());
							} else {
								writer.println(line);
							}
						}
						reader.close();
						writer.close();

						writer = new PrintWriter(new File(files[i], "pFind.bat"));
						writer.println("C:");
						writer.println("cd C:\\Users\\daz25hak\\MetaProAI\\Resources\\DBReducer\\pFind3\\bin");
						writer.println("searcher.exe " + (new File(files[i], "pFind.cfg")).getAbsolutePath());
						writer.close();
//					}
				} catch (IOException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
			}
		}
	}
	
	public static void getGenomeContribution(String in) {
		try {
			File[] files = (new File(in)).listFiles();
			for (int i = 0; i < files.length; i++) {
				HashMap<String, Integer> map = new HashMap<String, Integer>();
				File marker4HapFile = new File(files[i], "hap");
				if (marker4HapFile.exists()) {
					int count = 0;
					File selectFile = new File(marker4HapFile, "selected.tsv");
					BufferedReader reader = new BufferedReader(new FileReader(selectFile));
					String line = reader.readLine();
					while ((line = reader.readLine()) != null) {
						if (line.startsWith("#")) {
							break;
						} else {
							String[] cs = line.split("\t");
							if (cs.length == 2) {
								map.put(cs[0], count++);
							}
						}
					}
					reader.close();
				}
				File marker4File = new File(files[i], "marker4");
				HashMap<String, Integer> map2 = new HashMap<String, Integer>();
				if (marker4File.exists()) {
					File psmFile = new File(marker4File, "pFind-Filtered.spectra");
					if (psmFile.exists()) {
						PFindPsmReader reader = new PFindPsmReader(psmFile);
						PFindPSM[] psms = reader.getPsms();
						for (PFindPSM psm : psms) {
							String[] pros = psm.getProteins();
							if (pros.length == 1 && pros[0].startsWith("MGYG")) {
								String genome = pros[0].substring(0, pros[0].indexOf("_"));
								if (map2.containsKey(genome)) {
									map2.put(genome, map2.get(genome) + 1);
								} else {
									map2.put(genome, 1);
								}
							}
						}
					}
				}

				File testFile = new File(marker4File, "testgenome.tsv");
				PrintWriter writer = new PrintWriter(testFile);
				for (String genome : map.keySet()) {
					if (map2.containsKey(genome)) {
						writer.println(genome + "\t" + map.get(genome) + "\t" + map2.get(genome));
					}
				}
				writer.close();
			}

		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
	private static void compareTwoFile(String s1, String s2) {
		try {
			BufferedReader reader1 = new BufferedReader(new FileReader(s1));
			String line1 = null;

			BufferedReader reader2 = new BufferedReader(new FileReader(s2));
			String line2 = null;

			while (true) {
				line1 = reader1.readLine();
				line2 = reader2.readLine();
				if (line1 == null || line2 == null) {
					break;
				}
				
				if(!line1.equals(line2)) {
					System.out.println(line1);
					System.out.println(line2);
					System.out.println();
				}
			}

			reader1.close();
			reader2.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
	private static void getGenomeProteins(String in, String out, String eggnog) {

		try {
			HashMap<String, HashSet<String>> map = new HashMap<String, HashSet<String>>();
			BufferedReader reader = new BufferedReader(new FileReader(in));

			String line = reader.readLine();
			while ((line = reader.readLine()) != null) {
				String[] cs = line.split(",");
				String genome = cs[0].substring(0, cs[0].indexOf("_"));
				map.computeIfAbsent(genome, k -> new HashSet<>()).add(cs[0]);
			}
			reader.close();
			System.out.println(map.size());

			PrintWriter writer = new PrintWriter(out);
			boolean hasTitle = false;
			for (String genome : map.keySet()) {
				File file = new File(eggnog, genome + "_eggNOG.tsv");
				if (file.exists()) {
					HashSet<String> set = map.get(genome);
					System.out.println(genome+"\t"+set);
					
					reader = new BufferedReader(new FileReader(file));
					line = reader.readLine();
					if (!hasTitle) {
						writer.println(line);
						hasTitle = true;
					}
					while ((line = reader.readLine()) != null) {
						String[] cs = line.split("\t");
						if (set.contains(cs[0])) {
							writer.println(line);
						}
					}
					reader.close();
				}
			}
			writer.close();

		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
	private static void combineGenome(String in) {
		try {
			PrintWriter writer = new PrintWriter(new File(in, "CompMarker4SelectedGenome.tsv"));
			writer.println("File name\tRank\tGenome");
			File[] files = (new File(in)).listFiles();
			for (int i = 0; i < files.length; i++) {
				String name = files[i].getName();
				File hapFile = new File(files[i], "hap");
				File genomeFile = new File(hapFile, "comp_ranked.tsv");
				if (genomeFile.exists()) {
					int count = 0;
					BufferedReader reader = new BufferedReader(new FileReader(genomeFile));
					String line = reader.readLine();
					while ((line = reader.readLine()) != null) {
						String[] cs = line.split("\t");
//						if(cs[9].trim().length()>0) {
						if(cs[17].trim().length()>0) {
							writer.println(name + "\t" + (++count) + "\t" + cs[1]);
						} else {
							break;
						}
						
//						if (cs[4].equals("1.0")) {
//							break;
//						}
					}
					reader.close();
				}
			}
			writer.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
	private static void getProCount(String in) {
		try {
			int count = 0;
			File[] files = (new File(in)).listFiles();
			for(int i=0;i<files.length;i++) {
				if(files[i].getName().endsWith(".fasta")) {
					HashSet<String> genomeSet = new HashSet<String>();
					BufferedReader reader = new BufferedReader(new FileReader(files[i]));
					String line = reader.readLine();
					while ((line = reader.readLine()) != null) {
						if(line.startsWith(">")) {
							count++;
							int id = line.indexOf("_");
							if(id>1) {
								String genome = line.substring(1, line.indexOf("_"));
								genomeSet.add(genome);
							}else {
								System.out.println(line);
							}
							
						}
					}
					reader.close();
					System.out.println(files[i].getName()+"\t"+genomeSet.size());
				}
			}
			System.out.println(count);
			
		}catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
	private static void match300drugPeptide(String in, String out, String predicted) {
		try {
			BufferedReader reader = new BufferedReader(new FileReader(in));
			String line = reader.readLine();
			System.out.println(line);
			HashMap<String, HashSet<String>> peptideMap = new HashMap<String, HashSet<String>>();
			while ((line = reader.readLine()) != null) {
				String[] cs = line.split("\t");
				peptideMap.put(cs[0], new HashSet<String>());
			}
			reader.close();
			System.out.println(peptideMap.size());

			HashMap<String, Double> scoreMap = new HashMap<String, Double>();
			File[] files = (new File(predicted)).listFiles();
			for (File file : files) {
				reader = new BufferedReader(new FileReader(file));
				line = reader.readLine();
				while ((line = reader.readLine()) != null) {
					String[] cs = line.split("\t");
					if (peptideMap.containsKey(cs[1])) {
						peptideMap.get(cs[1]).add(cs[0]);
						double score = Double.parseDouble(cs[2]);
						if (scoreMap.containsKey(cs[1])) {
							if (score > scoreMap.get(cs[1])) {
								scoreMap.put(cs[1], score);
							}
						} else {
							scoreMap.put(cs[1], score);
						}
					}
				}
				reader.close();
			}

			PrintWriter writer = new PrintWriter(out);
			reader = new BufferedReader(new FileReader(in));

			String[] title = reader.readLine().split("\t");
			ArrayList<Integer> list = new ArrayList<Integer>();
			StringBuilder sb = new StringBuilder();
			sb.append("peptide\t");
			sb.append("proteins\t");
			sb.append("Detectability");
			for (int i = 0; i < title.length; i++) {
				if (title[i].startsWith("Intensity ")) {
					list.add(i);
					sb.append("\t").append(title[i]);
				}
			}
			writer.println(sb);

			while ((line = reader.readLine()) != null) {
				String[] cs = line.split("\t");
				if (peptideMap.containsKey(cs[0]) && peptideMap.get(cs[0]).size() > 0) {
					sb = new StringBuilder();
					sb.append(cs[0]).append("\t");
					for (String pro : peptideMap.get(cs[0])) {
						sb.append(pro).append(";");
					}
					sb.append("\t").append(scoreMap.get(cs[0]));
					for (int i = 0; i < list.size(); i++) {
						sb.append("\t").append(cs[list.get(i)]);
					}
					writer.println(sb);
				}
			}
			reader.close();
			writer.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
	private static void refine300Drug(String pepin, String genomein, String out) {
		try {
			HashSet<String> genomeSet = new HashSet<String>();
			BufferedReader reader = new BufferedReader(new FileReader(genomein));
			String line = reader.readLine();
			while ((line = reader.readLine()) != null) {
				if (line.startsWith("MGYG")) {
					String[] cs = line.split("\t");
					genomeSet.add(cs[0]);
				}
			}
			reader.close();
			System.out.println(genomeSet.size());

			PrintWriter writer = new PrintWriter(out);
			reader = new BufferedReader(new FileReader(pepin));
			line = reader.readLine();
			writer.println(line);
			while ((line = reader.readLine()) != null) {
				String[] cs = line.split("\t");
				String[] proteins = cs[1].split(";");
				ArrayList<String> prolist = new ArrayList<String>();
				for (String pro : proteins) {
					String genome = pro.substring(0, pro.indexOf("_"));
					if (genomeSet.contains(genome)) {
						prolist.add(pro);
					}
				}

				if (prolist.size() > 0) {
					StringBuilder sb = new StringBuilder();
					sb.append(cs[0]).append("\t");
					for (String pro : prolist) {
						sb.append(pro).append(";");
					}

					for (int i = 2; i < cs.length; i++) {
						sb.append("\t").append(cs[i]);
					}
					writer.println(sb);
				}
			}
			reader.close();
			writer.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
	private static void copyCfg(String in, String out) {
		File[] files = (new File(in)).listFiles();
		for (int i = 0; i < files.length; i++) {
			File hapFile = new File(files[i], "hap");
			if (hapFile.exists()) {
				File toFile = new File(out, files[i].getName());
				File toHapFile = new File(toFile, "hap");
				toHapFile.mkdirs();

				File[] hapFiles = hapFile.listFiles();
				for (int j = 0; j < hapFiles.length; j++) {
					String name = hapFiles[j].getName();
					if (name.startsWith("pFind")) {
						try {
							FileUtils.copyFile(hapFiles[j], new File(toHapFile, name));
						} catch (IOException e) {
							// TODO Auto-generated catch block
							e.printStackTrace();
						}
					}
				}
			}
		}
	}
	
	private static void sihumixTest(String in) {
		try {
			HashSet<String> set = new HashSet<String>();
			BufferedReader reader = new BufferedReader(new FileReader(in));
			String line = reader.readLine();
			String[] title = line.split("\t");
			int msms = -1;
			for (int i = 0; i < title.length; i++) {
				if (title[i].equals("MS/MS count S03")) {
//				if (title[i].equals("PS_SIHUMIx_S08")) {
					msms = i;
				}
			}
			while ((line = reader.readLine()) != null) {
				String[] cs = line.split("\t");
				if(!cs[msms].equals("0")) {
//					String seq = cs[1].replaceAll("I", "L");
					set.add(cs[1]);
				}
			}
			reader.close();

			System.out.println(set.size());
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

	public static void main(String[] args) {
		// TODO Auto-generated method stub
		PFindPsmReader.sihumixTest("C:\\Users\\daz25hak\\Documents\\data\\data\\marker4_focus.tsv");
//		PFindPsmReader.sihumixTest("C:\\Users\\daz25hak\\Documents\\data\\data\\SIHUMIx_peptide_paperdata.txt");
//		PFindPsmReader.copyCfg(
//				"Z:\\PRIDE_dataset\\IPX0003851000\\IPX0003851004\\DDA\\MetaPilot\\pfind_MAG_search_marker4_comp", 
//				"Z:\\PRIDE_dataset\\IPX0003851000\\IPX0003851004\\DDA\\MetaPilot\\pfind_MAG_search");
//		PFindPsmReader.match300drugPeptide("Z:\\PRIDE_dataset\\PXD056930_300drug\\search\\peptides.txt", 
//				"Z:\\PRIDE_dataset\\PXD056930_300drug\\search\\mgnify_peptides.txt", 
//				"C:\\Users\\daz25hak\\MetaProAI\\database\\MAG\\human-gut\\v2.0.2\\predicted");
//		PFindPsmReader.refine300Drug("Z:\\PRIDE_dataset\\PXD056930_300drug\\search\\mgnify_peptides.txt", 
//				"Z:\\PRIDE_dataset\\PXD056930_300drug\\selected.tsv", 
//				"Z:\\PRIDE_dataset\\PXD056930_300drug\\peptides.tsv");
//		PFindPsmReader.getProCount("Z:\\PRIDE_dataset\\PXD008738_Elo_DIA\\pool\\DIA\\MetaPilot\\marker_focus\\marker");
//		PFindPsmReader.cleanFolder("Z:\\Test\\BellAshley\\MetaPilot\\fragpipe_MAG_search");
//		PFindPsmReader.getProCount("C:\\Users\\daz25hak\\MetaProAI\\database\\MAG\\human-gut\\v2.0.2\\original_db");
		String folder = "Z:\\PRIDE_dataset\\PXD012724_rapidAim1\\MetaPilot\\";
		String file = "\\Leyuan_20170610_V1_ATRV\\hap\\pFind-Filtered.spectra";
//		PFindPsmReader.compare(folder + "pfind_MAG_search" + file, folder + "pfind_MAG_search_rib" + file);
		System.out.println();
		file = "\\Leyuan_20170610_V1_AZTH\\hap\\pFind-Filtered.spectra";
//		PFindPsmReader.compare(folder + "pfind_MAG_search" + file, folder + "pfind_MAG_search_rib" + file);
		System.out.println();
		file = "\\Leyuan_20170610_V1_BRBR\\hap\\pFind-Filtered.spectra";
//		PFindPsmReader.compare(folder + "pfind_MAG_search" + file, folder + "pfind_MAG_search_rib" + file);
//		PFindPsmReader.copyCfg("Z:\\PRIDE_dataset\\PXD023217_campi\\MetaPilot\\pfind_MAG_search_marker4_85");
//		PFindPsmReader.copyCfgFasta("Z:\\PRIDE_dataset\\PXD012724_rapidAim1\\MetaPilot\\pfind_MAG_search",
//				"C:\\Users\\daz25hak\\MetaProAI\\database\\MAG\\human-gut\\v2.0.2\\original_db");
//		PFindPsmReader.getGenomeContribution("Z:\\PRIDE_dataset\\PXD006118_3mix_type\\Run1\\MetaPilot\\pfind_MAG_search");
//		PFindPsmReader.compareFasta("Z:\\PRIDE_dataset\\PXD007819_xu_IBD\\MetaPilot\\pfind_MAG_search\\Xu_201606_HM406_TI\\hap"
//				+ "\\Xu_201606_HM406_TI.fasta", "Z:\\PRIDE_dataset\\PXD007819_xu_IBD\\sample_specific_db.fasta");
//		PFindPsmReader.compareFolder("Z:\\PRIDE_dataset\\PXD023217_campi\\MetaPilot\\pfind_MAG_search_marker2_comp", 
//				"Z:\\PRIDE_dataset\\PXD023217_campi\\MetaPilot\\pfind_MAG_search");
//		PFindPsmReader.compareTwoFile("Z:\\PRIDE_dataset\\PXD023217_campi\\SIHUMIx\\MetaPilot\\pfind_MAG_search\\S03\\marker4_70\\S03.fasta", 
//				"Z:\\PRIDE_dataset\\PXD023217_campi\\SIHUMIx\\MetaPilot\\pfind_MAG_search_rib\\S03\\hap\\S03.fasta");
//		PFindPsmReader.getGenomeProteins("C:\\Users\\daz25hak\\Documents\\data\\outputs\\Proteins_candidates.csv",
//				"C:\\Users\\daz25hak\\Documents\\data\\outputs\\Proteins_candidates_functions.tsv",
//				"C:\\Users\\daz25hak\\MetaProAI\\database\\MAG\\human-gut\\v2.0.2\\eggNOG");
//		PFindPsmReader.combineGenome("Z:\\PRIDE_dataset\\PXD023217_campi\\SIHUMIx\\MetaPilot\\pfind_MAG_search_focus");
		/*
		File[] files = (new File("Z:\\PRIDE_dataset\\PXD023217_campi\\SIHUMIx\\MetaPilot\\pfind_MAG_search")).listFiles();
		Arrays.sort(files);
		for (int i = 0; i < files.length; i++) {
			File hapFile = new File(files[i], "hap");
			if (hapFile.exists()) {
//				File hapResultFile = new File(hapFile, "pFind-Filtered.spectra");
				File hapResultFile = new File(files[i], "pFind-Filtered.spectra");
				if (hapResultFile.exists()) {
					PFindPsmReader reader = new PFindPsmReader(hapResultFile);
					HashSet<String> genomeSet = new HashSet<String>();
					for(PFindPSM psm: reader.getPsms()) {
						String[] proteins = psm.getProteins();
						for(String pro: proteins) {
							if(pro.startsWith("MGYG")) {
								String genome = pro.substring(0, pro.indexOf("_"));
								genomeSet.add(genome);
							}
						}
					}
					System.out.println(files[i].getName() + "\t" + reader.getPsms().length + "\t" + genomeSet.size());
				}
				
				File hapFastaFile = new File(new File(files[i], "hap"), files[i].getName() + ".fasta");
				try {
					int count = 0;
					HashSet<String> genomeSet = new HashSet<String>();
					BufferedReader reader = new BufferedReader(new FileReader(hapFastaFile));
					String line = null;
					while ((line = reader.readLine()) != null) {
						if (line.startsWith(">")) {
							count++;
							String genome = line.substring(1, line.indexOf("_"));
							genomeSet.add(genome);
						}
					}
					reader.close();
					System.out.println(files[i].getName() + "\t" + count + "\t" + genomeSet.size());
				} catch (IOException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
				
				hapFastaFile = new File(new File(files[i], "marker4"), files[i].getName() + ".fasta");
				try {
					int count = 0;
					HashSet<String> genomeSet = new HashSet<String>();
					BufferedReader reader = new BufferedReader(new FileReader(hapFastaFile));
					String line = null;
					while ((line = reader.readLine()) != null) {
						if (line.startsWith(">")) {
							count++;
							String genome = line.substring(1, line.indexOf("_"));
							genomeSet.add(genome);
						}
					}
					reader.close();
					System.out.println(files[i].getName() + "\t" + count + "\t" + genomeSet.size());
				} catch (IOException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
				
				hapFastaFile = new File(new File(files[i], "marker4_70"), files[i].getName() + ".fasta");
				try {
					int count = 0;
					HashSet<String> genomeSet = new HashSet<String>();
					BufferedReader reader = new BufferedReader(new FileReader(hapFastaFile));
					String line = null;
					while ((line = reader.readLine()) != null) {
						if (line.startsWith(">")) {
							count++;
							String genome = line.substring(1, line.indexOf("_"));
							genomeSet.add(genome);
						}
					}
					reader.close();
					System.out.println(files[i].getName() + "\t" + count + "\t" + genomeSet.size());
				} catch (IOException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
				
			}
		}
		*/
		/*
		PFindPsmReader reader = new PFindPsmReader("Z:\\PRIDE_dataset\\PXD023217_campi\\SIHUMIx\\MetaPilot\\pfind_MAG_search_old\\S08"
				+ "\\pFind-Filtered.spectra");
		System.out.println(reader.getPsms().length);
		HashSet<String> set = new HashSet<String>();
		for(PFindPSM psm: reader.getPsms()) {
			String[] proteins = psm.getProteins();
			for(String pro: proteins) {
				String genome = pro.substring(0, pro.indexOf("_"));
				set.add(genome);
				
			}
		}
		System.out.println(set);
		*/
	}
}
