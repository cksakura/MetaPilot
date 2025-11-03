package uk.ac.quadram.metapilot.dbSearch.dbReducer;

import java.io.File;
import java.io.IOException;

import uk.ac.quadram.metapilot.dbSearch.pfind.PFindConfig;

public class DBReducerConfig extends PFindConfig {

	public DBReducerConfig(File cfgFile) {
		super(cfgFile);
		// TODO Auto-generated constructor stub
	}

	protected void parse(File cfgFile) {
		try {
			findMap = this.parseMap(cfgFile);
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

	public File searchMgf(String[] raws, String result, String db, String ms2Type, double psmFdr, String[] fixedMods,
			String[] variMods) {

		StringBuilder fixsb = new StringBuilder();
		for (int i = 0; i < fixedMods.length; i++) {
			String mod = fixedMods[i];
			int spaceIndex = mod.indexOf(" ");
			if (spaceIndex > 0) {
				fixsb.append(mod.substring(0, spaceIndex)).append(";");
			} else {
				// If no space found, use the entire string
				fixsb.append(mod).append(";");
			}
		}

		StringBuilder varisb = new StringBuilder();
		for (int i = 0; i < variMods.length; i++) {
			String mod = variMods[i];
			int spaceIndex = mod.indexOf(" ");
			if (spaceIndex > 0) {
				varisb.append(mod.substring(0, spaceIndex)).append(";");
			} else {
				// If no space found, use the entire string
				varisb.append(mod).append(";");
			}
		}
		if (varisb.length() > 0) {
			varisb.deleteCharAt(varisb.length() - 1);
		}

		this.findMap.get("fixmod").setValue(fixsb.toString());
		this.findMap.get("selectmod").setValue(varisb.toString());
		this.findMap.get("ProteinDatabase").setValue(db);
		this.findMap.get("SearchRes").setValue(result);
		this.findMap.get("threshold").setValue(String.valueOf(psmFdr));

		this.findMap.get("activation_type").setValue(ms2Type);
		this.findMap.get("msmsnum").setValue(String.valueOf(raws.length));

		this.findMap.get("msmstype").setValue("MGF");

		StringBuilder findRawSb = new StringBuilder(raws[0]);

		for (int i = 1; i < raws.length; i++) {
			findRawSb.append("\n");
			findRawSb.append("msmspath").append(i + 1).append("=").append(raws[i]);
		}

		this.findMap.get("msmspath1").setValue(findRawSb.toString());

		try {
			this.pFindCfg = new File(result, "DBReducer.cfg");
			this.export(pFindCfg, findMap);
		} catch (IOException e) {
			// TODO Auto-generated method stub
			e.printStackTrace();
		}

		return pFindCfg;
	}

}
