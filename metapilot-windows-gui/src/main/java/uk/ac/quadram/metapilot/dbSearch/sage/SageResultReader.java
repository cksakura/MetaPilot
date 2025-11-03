package uk.ac.quadram.metapilot.dbSearch.sage;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;

import uk.ac.quadram.metapilot.core.math.MathTool;
import uk.ac.quadram.metapilot.task.io.pep.AbstractMetaPeptideReader;

public class SageResultReader extends AbstractMetaPeptideReader {

	private SagePeptide[] peptides;
	private SagePSM[] psms;
	private String[] intensityTitles;

	public SageResultReader(String file) {
		this(new File(file));
	}

	public SageResultReader(File file) {
		super(file);
		read();
	}

	private void read() {
		File psmFile = new File(this.getFile(), "results.sage.tsv");
		SagePsmReader reader = new SagePsmReader(psmFile);
		HashMap<String, ArrayList<SagePSM>> psmMap = new HashMap<String, ArrayList<SagePSM>>();
		for (SagePSM psm : reader.getPsms()) {
			String modseq = psm.getModSeq();
			if (psmMap.containsKey(modseq)) {
				psmMap.get(modseq).add(psm);
			} else {
				ArrayList<SagePSM> psmList = new ArrayList<SagePSM>();
				psmList.add(psm);
				psmMap.put(modseq, psmList);
			}
		}

		HashSet<String> modPepSet = new HashSet<String>();
		ArrayList<SagePeptide> peptideList = new ArrayList<SagePeptide>();
		try (BufferedReader br = new BufferedReader(new FileReader(new File(this.getFile(), "lfq.tsv")))) {
			ArrayList<String> list = new ArrayList<String>();
			String line = br.readLine();
			String[] title = line.split("\t");
			int pepId = -1;
			int proId = -1;
			int qvalueId = -1;
			int scoreId = -1;
			ArrayList<Integer> idList = new ArrayList<Integer>();

			for (int i = 0; i < title.length; i++) {
				if (title[i].equals("peptide")) {
					pepId = i;
				} else if (title[i].equals("proteins")) {
					proId = i;
				} else if (title[i].equals("q_value")) {
					qvalueId = i;
				} else if (title[i].equals("score")) {
					scoreId = i;
				} else {
					if (title[i].endsWith(".mzML")) {
						String fileName = title[i].substring(0, title[i].length() - ".mzML".length());
						list.add(fileName);
						idList.add(i);
					}
				}
			}
			this.intensityTitles = list.toArray(String[]::new);
			int[] ids = new int[idList.size()];
			for (int i = 0; i < ids.length; i++) {
				ids[i] = idList.get(i);
			}

			while ((line = br.readLine()) != null) {
				String[] cs = line.split("\t");
				String modseq = cs[pepId];
				if (psmMap.containsKey(modseq)) {
					modPepSet.add(modseq);
					ArrayList<SagePSM> psmList = psmMap.get(modseq);
					SagePSM psm0 = psmList.get(0);
					String sequence = psm0.getSequence();
					String[] proteins = cs[proId].split(";");
					double[] intensity = new double[ids.length];
					int[] idenType = new int[ids.length];
					int[] ms2Count = new int[ids.length];
					for (int i = 0; i < intensity.length; i++) {
						intensity[i] = Double.parseDouble(cs[ids[i]]);
						if (intensity[i] > 0) {
							idenType[i] = 1;
						} else {
							idenType[i] = 2;
						}
					}

					double[] indiPEP = new double[ids.length];
					Arrays.fill(indiPEP, 1);
					for (SagePSM psm : psmList) {
						String fileName = psm.getFileName();
						for (int i = 0; i < intensityTitles.length; i++) {
							if (fileName.equals(intensityTitles[i])) {
								idenType[i] = 0;
								ms2Count[i]++;
								if (indiPEP[i] < psm.getqValue()) {
									indiPEP[i] = psm.getqValue();
								}
							}
						}
					}
					SagePeptide peptide = new SagePeptide(sequence, modseq, proteins, intensity, idenType,
							Double.parseDouble(cs[qvalueId]), Double.parseDouble(cs[scoreId]));
					HashSet<Integer> chargeSet = new HashSet<Integer>();

					for (int i = 0; i < psmList.size(); i++) {
						SagePSM psm = psmList.get(i);
						peptide.addPSM(psm);
						chargeSet.add(psm.getCharge());
					}
					int[] charges = new int[chargeSet.size()];
					int chargeId = 0;
					for (Integer charge : chargeSet) {
						charges[chargeId++] = charge;
					}
					peptide.setCharges(charges);
					peptide.setMs2Counts(ms2Count);
					peptide.setIndividualPEP(indiPEP);
					peptideList.add(peptide);
				}
			}
			br.close();

		} catch (IOException e) {

		}

		this.peptides = peptideList.toArray(SagePeptide[]::new);

		ArrayList<SagePSM> psmList = new ArrayList<SagePSM>();
		for (SagePSM psm : reader.getPsms()) {
			String modseq = psm.getModSeq();
			if (modPepSet.contains(modseq)) {
				psmList.add(psm);
			}
		}

		this.psms = psmList.toArray(SagePSM[]::new);
	}

	@Override
	public String getQuanMode() {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public SagePeptide[] getMetaPeptides() {
		// TODO Auto-generated method stub
		return peptides;
	}

	public SagePSM[] getPSMs() {
		return psms;
	}

	@Override
	public Object[] getTitleObjs() {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public String[] getIntensityTitle() {
		// TODO Auto-generated method stub
		return intensityTitles;
	}
	
	public static void main(String[] args) {
		SageResultReader reader = new SageResultReader(
				"Z:\\PRIDE_dataset\\PXD050781\\DDA21\\MetaPilot\\sage_MAG_search\\combined");
		SagePeptide[] peptides = reader.getMetaPeptides();
		for(int i=0;i<peptides.length;i++) {
			double[] pepIntensity = peptides[i].getIntensity();
			double totalIntensity = MathTool.getTotal(pepIntensity);
			int[] idenType = peptides[i].getIdenType();
			int idenCount = 0;
			for (int j = 0; j < idenType.length; j++) {
				if (idenType[j] == 0) {
					idenCount++;
				}
			}
			System.out.println(totalIntensity+"\t"+idenCount);
		}
	}

}
