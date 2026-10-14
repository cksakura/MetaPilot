package uk.ac.quadram.metapilot.quant.flashLFQ;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Date;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class FlashLfqPsmReader {

	private FlashLfqPsm[] psms;
	
	private SimpleDateFormat format = new SimpleDateFormat("MM_dd_yyyy_HH_mm_ss");
	private static final Logger LOGGER = LogManager.getLogger(FlashLfqPsmReader.class);

	public FlashLfqPsmReader(String file) throws IOException {
		this(new File(file));
	}

	public FlashLfqPsmReader(File file) throws IOException {
		this.read(file);
	}

	private void read(File file) {
		ArrayList<FlashLfqPsm> list = new ArrayList<FlashLfqPsm>();
		try (BufferedReader reader = new BufferedReader(new FileReader(file))) {
			String line = reader.readLine();
			String[] title = line.split("\t");

			int fileNameId = -1;
			int scanId = -1;
			int baseSequenceId = -1;
			int fullSequenceId = -1;
			int expMassId = -1;
			int rtId = -1;
			int chargeId = -1;
			int proteinsId = -1;
			int scoreId = -1;
			int PEPId = -1;
			int calMassId = -1;
			int massDiffId = -1;
			int missId = -1;

			for (int i = 0; i < title.length; i++) {
				if (title[i].equals("File Name")) {
					fileNameId = i;
				} else if (title[i].equals("Scan number")) {
					scanId = i;
				} else if (title[i].equals("Base Sequence")) {
					baseSequenceId = i;
				} else if (title[i].equals("Full Sequence")) {
					fullSequenceId = i;
				} else if (title[i].equals("Peptide Monoisotopic Mass")) {
					// FlashLFQ's generic format: the THEORETICAL mass it searches MS1 for.
					// OpenPsmTsvWriter.addPFindPSM writes calMH - proton here and the
					// experimental mass under "Peptide Mass"; this reader used to map them
					// the other way round, so getCalMass() returned the experimental mass
					// (12mix pFind psms.tsv: 2918.524881 here vs 2918.528629 under "Peptide
					// Mass", "Mass difference" 0.003748 = the latter minus the former).
					calMassId = i;
				} else if (title[i].equals("Scan Retention Time")) {
					rtId = i;
				} else if (title[i].equals("Precursor Charge")) {
					chargeId = i;
				} else if (title[i].equals("Protein Accession")) {
					proteinsId = i;
				} else if (title[i].equals("Score")) {
					scoreId = i;
				} else if (title[i].equals("PEP")) {
					PEPId = i;
				} else if (title[i].equals("Peptide Mass")) {
					expMassId = i;
				} else if (title[i].equals("Mass difference")) {
					massDiffId = i;
				} else if (title[i].equals("Missed cleavages")) {
					missId = i;
				}
			}

			while ((line = reader.readLine()) != null) {
				String[] cs = line.split("\t");
				FlashLfqPsm psm = new FlashLfqPsm(cs[fileNameId], Integer.parseInt(cs[scanId]), cs[baseSequenceId], cs[fullSequenceId],
						Double.parseDouble(cs[expMassId]), Double.parseDouble(cs[rtId]), Integer.parseInt(cs[chargeId]),
						cs[proteinsId].split(";"), Double.parseDouble(cs[scoreId]), Double.parseDouble(cs[PEPId]),
						Double.parseDouble(cs[calMassId]), Double.parseDouble(cs[massDiffId]),
						Integer.parseInt(cs[missId]));
				list.add(psm);
			}

			reader.close();
		} catch (IOException e) {
			System.out.println(format.format(new Date()) + "\t" + ": error in reading PSMs from " + file);
			LOGGER.error("Error in reading PSMs from " + file, e);
		}

		this.psms = list.toArray(new FlashLfqPsm[list.size()]);
	}

	public FlashLfqPsm[] getPSMs() {
		return psms;
	}

	public static void main(String[] args) {
		// TODO Auto-generated method stub

	}

}
