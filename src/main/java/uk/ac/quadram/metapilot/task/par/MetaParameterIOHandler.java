package uk.ac.quadram.metapilot.task.par;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.json.JSONArray;
import org.json.JSONException;
import org.json.JSONObject;

import uk.ac.quadram.metapilot.core.mod.IsobaricTag;
import uk.ac.quadram.metapilot.core.taxonomy.RootType;
import uk.ac.quadram.metapilot.core.taxonomy.TaxonomyRanks;
import uk.ac.quadram.metapilot.dbSearch.MaxQuant.MaxquantModification;
import uk.ac.quadram.metapilot.task.MetaLabWorkflowType;
import uk.ac.quadram.metapilot.task.pfind.par.MetaParaIOPFind;
import uk.ac.quadram.metapilot.task.pfind.par.MetaParameterPFind;
import uk.ac.quadram.metapilot.task.v2.par.MaxQuantModIO;
import uk.ac.quadram.metapilot.task.v2.par.MetaOpenPar;
import uk.ac.quadram.metapilot.task.v2.par.MetaParaIOMQ;

public class MetaParameterIOHandler {

	private static Logger LOGGER = LogManager.getLogger(MetaParameterIOHandler.class);
	
	public static MetaLabWorkflowType getWorkflow(String json) {
		return getWorkflow(new File(json));
	}
	
	public static MetaLabWorkflowType getWorkflow(File json) {
		if (!json.exists() || json.length() == 0) {
			LOGGER.error("Error in reading MetaLab parameter file " + json);
			return null;
		}

		StringBuilder sb = new StringBuilder();
		try (BufferedReader reader = new BufferedReader(new FileReader(json))) {
			String line = null;
			while ((line = reader.readLine()) != null) {
				sb.append(line);
			}
			reader.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			LOGGER.error("Error in reading MetaLab parameter file " + json, e);
			return null;
		}

		JSONObject obj = null;
		try {
			obj = new JSONObject(sb.toString());
		} catch (JSONException e) {
			LOGGER.error("Error in reading MetaLab parameter file " + json, e);
			return null;
		}

		if (obj.has("workflowType")) {
			MetaLabWorkflowType workflow = MetaLabWorkflowType.valueOf(obj.getString("workflowType"));
			return workflow;
		} else {
			LOGGER.error("There is no workflowType information in " + json);
			return null;
		}
	}

	public static MetaBaseParameter parse(String json) {
		return parse(new File(json));
	}

	public static MetaBaseParameter parse(File json) {

		if (!json.exists() || json.length() == 0) {
			MetaParaIOMQ.exportBlank(json);
		}

		StringBuilder sb = new StringBuilder();
		try (BufferedReader reader = new BufferedReader(new FileReader(json))) {
			String line = null;
			while ((line = reader.readLine()) != null) {
				sb.append(line);
			}
			reader.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			LOGGER.error("Error in reading MetaLab parameter file " + json, e);
			return null;
		}

		JSONObject obj = null;
		try {
			obj = new JSONObject(sb.toString());
		} catch (JSONException e) {
			LOGGER.error("Error in reading MetaLab parameter file " + json, e);
			return null;
		}

		return parse(obj);
	}

	public static MetaBaseParameter parseMain(JSONObject obj) {
		String workflowType = obj.has("workflowType") ? obj.getString("workflowType") : "MaxQuantWorkflow";
		JSONArray rawsExpArray = obj.getJSONArray("rawExpName");
		String[] rawFiles = new String[rawsExpArray.length()];
		String[] expNames = new String[rawsExpArray.length()];
		int[] fraction = new int[rawsExpArray.length()];
		int[] replicate = new int[rawsExpArray.length()];

		for (int i = 0; i < rawFiles.length; i++) {
			JSONObject obji = rawsExpArray.getJSONObject(i);
			rawFiles[i] = obji.getString("path");
			expNames[i] = obji.getString("experiment");
			if (obji.has("fraction")) {
				fraction[i] = obji.getInt("fraction");
			}
			if (obji.has("replicate")) {
				replicate[i] = obji.getInt("replicate");
			}
		}

		MetaData metadata;

		int metaTypeCount = obj.getInt("metaCount");
		
		String[] labelTitle = new String[] {};
		if(obj.has("labelTitle")) {
			JSONArray labelTitleArray = obj.getJSONArray("labelTitle");
			labelTitle = new String[labelTitleArray.length()];
			for (int i = 0; i < labelTitle.length; i++) {
				labelTitle[i] = labelTitleArray.getJSONObject(i).getString("label");
			}
		}
		
		int labelCount = labelTitle.length;

		String[][] metaInfo;
		if (labelCount == 0) {
			if (metaTypeCount == 0) {
				metadata = new MetaData(rawFiles, expNames, fraction, replicate);
			} else {
				metaInfo = new String[rawFiles.length][metaTypeCount];
				JSONArray metainfoArray = obj.getJSONArray("metainfo");
				for (int i = 0; i < metaInfo.length; i++) {
					JSONObject obji = metainfoArray.getJSONObject(i);
					for (int j = 0; j < metaTypeCount; j++) {
						metaInfo[i][j] = obji.getString("meta " + (j + 1));
					}
				}
				metadata = new MetaData(rawFiles, expNames, fraction, replicate, metaTypeCount, metaInfo);
			}

		} else {
			String[] labelNames = new String[expNames.length * labelCount];
			JSONArray labelArray = obj.getJSONArray("labelExpName");
			for (int i = 0; i < expNames.length; i++) {
				for (int j = 0; j < labelCount; j++) {
					int id = i * labelCount + j;
					JSONObject obji = labelArray.getJSONObject(id);
					labelNames[id] = obji.getString(expNames[i] + " " + labelTitle[j]);
				}
			}
			if (metaTypeCount == 0) {
				metadata = new MetaData(rawFiles, expNames, fraction, replicate, metaTypeCount, labelTitle, labelNames);
			} else {
				metaInfo = new String[rawFiles.length * labelCount][metaTypeCount];
				JSONArray metainfoArray = obj.getJSONArray("metainfo");
				for (int i = 0; i < metaInfo.length; i++) {
					JSONObject obji = metainfoArray.getJSONObject(i);
					for (int j = 0; j < metaTypeCount; j++) {
						metaInfo[i][j] = obji.getString("meta " + (j + 1));
					}
				}
				metadata = new MetaData(rawFiles, expNames, fraction, replicate, metaTypeCount, metaInfo, labelTitle,
						labelNames);
			}

			JSONArray isobaricRefArray = obj.getJSONArray("isobaricRefs");
			String[] isobaricRefs = new String[isobaricRefArray.length()];
			for (int i = 0; i < isobaricRefArray.length(); i++) {
				isobaricRefs[i] = isobaricRefArray.getJSONObject(i).getString("refName");
			}

			JSONArray selectRefArray = obj.getJSONArray("selectRefs");
			boolean[] selectRefs = new boolean[selectRefArray.length()];
			for (int i = 0; i < selectRefArray.length(); i++) {
				selectRefs[i] = selectRefArray.getJSONObject(i).getBoolean("select");
			}

			metadata.setIsobaricReference(isobaricRefs);
			metadata.setSelectRef(selectRefs);

			int refChannelId = obj.has("refChannelId") ? obj.optInt("refChannelId") : -1;
			metadata.setRefChannelId(refChannelId);
		}

		String result = obj.getString("result");

		String microDb = obj.getString("microDb");

		boolean appendHostDb = obj.getBoolean("appendHostDb");

		String hostDb = obj.getString("hostDb");

		String ms2ScanMode = obj.has("MS2ScanMode")?  obj.getString("MS2ScanMode"): "";

		int coreCount = obj.optInt("coreCount") == 0 ? 1 : obj.optInt("coreCount");

		int threadCount = obj.getInt("threadCount");

		boolean isMetaWorkflow = obj.getBoolean("isMetaWorkflow");

		MetaSsdbCreatePar mscp;
		if (obj.has("ssdb")) {
			boolean ssdb = obj.getBoolean("ssdb");
			boolean ssopen = obj.getBoolean("ssOpen");
			boolean cluster = obj.getBoolean("cluster");
			double xtandemFDR = obj.getDouble("xtandemFDR");
			double xtandemEvalue = obj.getDouble("xtandemEvalue");
			mscp = new MetaSsdbCreatePar(ssdb, cluster, ssopen, xtandemFDR, xtandemEvalue);
		} else {
			mscp = MetaSsdbCreatePar.getDefault();
		}

		MetaBaseParameter metaPar = new MetaBaseParameter(workflowType, metadata, result, microDb, appendHostDb, hostDb,
				ms2ScanMode, coreCount, threadCount, isMetaWorkflow, mscp);

		return metaPar;
	}

	public static MetaBaseParameter parse(JSONObject obj) {

		MetaBaseParameter metaPar;

		String workflowType = obj.has("workflowType") ? obj.getString("workflowType") : "MaxQuantWorkflow";

		MetaLabWorkflowType mwf = MetaLabWorkflowType.valueOf(workflowType);
		if (mwf == MetaLabWorkflowType.pFindWorkflow) {
			metaPar = MetaParaIOPFind.parse(obj);
		} else {
			metaPar = MetaParaIOMQ.parse(obj);
		}

		return metaPar;
	}
	
	public <T extends MetaBaseParameter> T parse(File file, Class<T> parameterType) {
		if (!file.exists() || file.length() == 0) {
//			exportBlank(file);
		}

		StringBuilder sb = new StringBuilder();
		String line = null;
		try (BufferedReader reader = new BufferedReader(new FileReader(file))) {
			while ((line = reader.readLine()) != null) {
				sb.append(line);
			}
			reader.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			LOGGER.error("Error in reading MetaPilot parameter file " + file, e);
			return null;
		}

		JSONObject obj = null;
		try {
			obj = new JSONObject(sb.toString());
		} catch (JSONException e) {
			LOGGER.error("Error in reading MetaPilot parameter file " + file, e);
			return null;
		}
		MetaBaseParameter metaPar = MetaParameterIOHandler.parseMain(obj);

		if (parameterType == MetaParameterPFind.class) {
			ArrayList<String> fixModList = new ArrayList<String>();
			JSONArray fixArray = obj.getJSONArray("fixMods");
			if (fixArray != null) {
				for (int i = 0; i < fixArray.length(); i++) {
					String name = fixArray.getJSONObject(i).getString("name");
					fixModList.add(name);
				}
			}
			String[] fixMods = fixModList.toArray(new String[fixModList.size()]);

			ArrayList<String> variModList = new ArrayList<String>();
			JSONArray variArray = obj.getJSONArray("variMods");
			if (variArray != null) {
				for (int i = 0; i < variArray.length(); i++) {
					String name = variArray.getJSONObject(i).getString("name");
					variModList.add(name);
				}
			}
			String[] variMods = variModList.toArray(new String[variModList.size()]);

			String enzymeString = obj.has("enzyme") ? obj.getString("enzyme") : "Trypsin";

			int missCleavages = obj.has("missCleavages") ? obj.getInt("missCleavages") : 2;

			int digestMode = obj.has("digestMode") ? obj.getInt("digestMode") : 0;

			String quanMode = obj.has("quanMode") ? obj.getString("quanMode") : MetaConstants.labelFree;

			boolean combineLabel = obj.has("combineLabel") ? obj.getBoolean("combineLabel") : false;

			HashMap<String, MaxquantModification> isobaricMap = MaxQuantModIO.getIsobaricMap();

			ArrayList<MaxquantModification> isobaricList = new ArrayList<MaxquantModification>();

			IsobaricTag isobaricTag = null;
			if (obj.has("isobaricTag")) {
				String tag = obj.getString("isobaricTag");
				if (tag.length() > 0) {
					isobaricTag = IsobaricTag.valueOf(tag);
				}
			}

			double[][] isoCorFactor = new double[][] {};
			JSONArray isobaricArray = obj.getJSONArray("isobaric");
			if (isobaricArray != null) {
				isoCorFactor = new double[isobaricArray.length()][4];
				for (int i = 0; i < isobaricArray.length(); i++) {
					String[] values = isobaricArray.getJSONObject(i).getString("name").split("_");
					if (isobaricMap.containsKey(values[0])) {
						isobaricList.add(isobaricMap.get(values[0]));
						for (int j = 0; j < 4; j++) {
							isoCorFactor[i][j] = Double.parseDouble(values[j + 1]);
						}
					}
				}
			}
			MaxquantModification[] isobaric = isobaricList.toArray(new MaxquantModification[isobaricList.size()]);

			boolean openSearch = obj.has("isOpenSearch") ? obj.getBoolean("isOpenSearch") : true;
			double openPsmFDR = obj.getDouble("openPsmFDR");
			double openProFDR = obj.getDouble("openProFDR");
			boolean matchBetweenRuns = obj.has("mbr") ? obj.getBoolean("mbr") : true;

			MetaOpenPar mosp = new MetaOpenPar(openSearch, openPsmFDR, openProFDR, matchBetweenRuns);

			boolean buildIn = obj.getBoolean("buildIn");
			boolean unipept = obj.getBoolean("unipept");

			JSONArray usedRootArray = obj.getJSONArray("usedRoots");
			HashSet<RootType> usedRootSet = new HashSet<RootType>();
			for (int i = 0; i < usedRootArray.length(); i++) {
				String name = usedRootArray.getJSONObject(i).getString("root");
				usedRootSet.add(RootType.valueOf(name));
			}

			TaxonomyRanks ignoreBlankRank = TaxonomyRanks.valueOf(obj.getString("ignoreBlankRank"));

			JSONArray ignoreTaxaArray = obj.getJSONArray("excludeTaxa");
			HashSet<String> ignoreTaxaSet = new HashSet<String>();
			for (int i = 0; i < ignoreTaxaArray.length(); i++) {
				ignoreTaxaSet.add(ignoreTaxaArray.getJSONObject(i).getString("excludeTaxon"));
			}
			int leastPepCount = obj.getInt("leastPepCount");

			MetaPep2TaxaPar mptp = new MetaPep2TaxaPar(buildIn, unipept, usedRootSet, ignoreBlankRank, ignoreTaxaSet,
					leastPepCount);

//			MetaParameterPFind par = new MetaParameterPFind(metaPar, fixMods, variMods, enzymeString, missCleavages,
//					digestMode, quanMode, isobaric, isobaricTag, combineLabel, isoCorFactor, openSearch, mosp, mptp);
//			return par;
		}
		return null;
	}
	
	 // Type-specific convenience methods
    public MetaParameterPFind parsePFindParameter(File file) {
        return parse(file, MetaParameterPFind.class);
    }
}
