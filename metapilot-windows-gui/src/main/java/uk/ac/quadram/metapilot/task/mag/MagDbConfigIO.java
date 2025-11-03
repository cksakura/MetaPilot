package uk.ac.quadram.metapilot.task.mag;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.text.ParseException;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Date;
import java.util.HashMap;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.json.JSONArray;
import org.json.JSONException;
import org.json.JSONObject;
import org.json.JSONWriter;

import uk.ac.quadram.metapilot.core.tools.FormatTool;
import uk.ac.quadram.metapilot.task.mag.MagDbItem.MagModel;
import uk.ac.quadram.metapilot.task.mag.MagDbItem.MagPepLib;
import uk.ac.quadram.metapilot.task.mag.MagDbItem.MagVersion;

public class MagDbConfigIO {

	public static final String version = "MAGDB 1.0";
	public static final String versionFile = "MAGDB_1_0.json";

	private static String MGnify = "MGnify";
	private static String localMAGs = "localMAGs";
	private static SimpleDateFormat sdf = FormatTool.getDateFormat();
	private static Logger LOGGER = LogManager.getLogger(MagDbConfigIO.class);

	public static MagDbItem[][] parse(String json) {
		return parse(new File(json));
	}

	public static MagDbItem[][] parse(File json) {
		if (!json.exists() || json.length() == 0) {
			exportBlank(json);
		}

		BufferedReader reader = null;
		try {
			reader = new BufferedReader(new FileReader(json));
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			LOGGER.error("Error in reading MetaLab parameter file " + json, e);
			return null;
		}

		StringBuilder sb = new StringBuilder();
		String line = null;
		try {
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

	public static MagDbItem[][] parse(JSONObject obj) {

		JSONArray mgnifyArray = obj.getJSONArray(MGnify);
		ArrayList<MagDbItem> mgnifyList = new ArrayList<MagDbItem>();
		if (mgnifyArray != null) {
			for (int i = 0; i < mgnifyArray.length(); i++) {
				JSONObject obji = mgnifyArray.getJSONObject(i);

				String catalogueID = obji.getString("catalogueID");
				String identifier = obji.getString("identifier");
				File repoFile = obji.isNull("repository") ? null : new File(obji.getString("repository"));

				MagDbItem magDbItem = new MagDbItem(catalogueID, identifier, repoFile);

				JSONArray versionArray = obji.getJSONArray("versionlist");
				for (int j = 0; j < versionArray.length(); j++) {
					JSONObject objj = versionArray.getJSONObject(j);

					String catalogueVersion = objj.getString("version");
					int speciesCount = objj.getInt("speciesCount");
					Date date = null;
					try {
						date = sdf.parse(objj.getString("date"));
					} catch (JSONException | ParseException e) {
						// TODO Auto-generated catch block
						e.printStackTrace();
					}
					boolean localAvailable = objj.getBoolean("localAvailable");
					String dbName = objj.getString("dbName");
					String funcName = objj.getString("funcName");
					String taxaName = objj.getString("taxaName");
					magDbItem.addVersion(catalogueVersion, speciesCount, date, localAvailable, dbName, funcName,
							taxaName);
				}

				JSONArray mgygArray = obji.getJSONArray("MGYG");
				for (int j = 0; j < mgygArray.length(); j++) {
					JSONObject objj = mgygArray.getJSONObject(j);
					String parent = objj.getString("parent");

					JSONArray childrenArray = objj.getJSONArray("children");
					for (int k = 0; k < childrenArray.length(); k++) {
						JSONObject objk = childrenArray.getJSONObject(k);
						String child = objk.getString("child");
						magDbItem.addMgygName(parent, child);
					}
				}

				if (obji.has("models")) {
					JSONArray modelArray = obji.getJSONArray("models");
					for (int j = 0; j < modelArray.length(); j++) {
						JSONObject objj = modelArray.getJSONObject(j);
						String modelName = objj.getString("modelName");
						int trainSize = objj.getInt("trainSize");
						int modelSize = objj.getInt("modelSize");
						double corr = objj.getDouble("corr");
						double med = objj.getDouble("med");
						magDbItem.addModel(modelName, trainSize, modelSize, corr, med);
					}
				}

				if (obji.has("pepLibs")) {
					JSONArray pepLibArray = obji.getJSONArray("pepLibs");
					for (int j = 0; j < pepLibArray.length(); j++) {
						JSONObject objj = pepLibArray.getJSONObject(j);
						String pepLibName = objj.getString("pepLibName");
						int pepCount = objj.getInt("pepCount");
						String protease = objj.getString("protease");
						int miss = objj.getInt("miss");
						magDbItem.addPepLib(pepLibName, pepCount, protease, miss);
					}
				}

				mgnifyList.add(magDbItem);
			}
		}

		JSONArray localArray = obj.getJSONArray(localMAGs);
		ArrayList<MagDbItem> localList = new ArrayList<MagDbItem>();
		if (localArray != null) {
			for (int i = 0; i < localArray.length(); i++) {
				JSONObject obji = localArray.getJSONObject(i);

				String catalogueID = obji.getString("catalogueID");
				String identifier = obji.getString("identifier");
				File repoFile = obji.isNull("repository") ? null : new File(obji.getString("repository"));
				MagDbItem magDbItem = new MagDbItem(catalogueID, identifier, repoFile);

				JSONArray versionArray = obji.getJSONArray("versionlist");
				for (int j = 0; j < versionArray.length(); j++) {
					JSONObject objj = versionArray.getJSONObject(j);

					String catalogueVersion = objj.getString("version");
					int speciesCount = objj.getInt("speciesCount");
					Date date = null;
					try {
						date = sdf.parse(objj.getString("date"));
					} catch (JSONException | ParseException e) {
						// TODO Auto-generated catch block
						e.printStackTrace();
					}
					boolean localAvailable = objj.getBoolean("localAvailable");
					String dbName = objj.getString("dbName");
					String funcName = objj.getString("funcName");
					String taxaName = objj.getString("taxaName");
					magDbItem.addVersion(catalogueVersion, speciesCount, date, localAvailable, dbName, funcName,
							taxaName);
				}

				if (obji.has("models")) {
					JSONArray modelArray = obji.getJSONArray("models");
					for (int j = 0; j < modelArray.length(); j++) {
						JSONObject objj = modelArray.getJSONObject(j);
						String modelName = objj.getString("modelName");
						int trainSize = objj.getInt("trainSize");
						int modelSize = objj.getInt("modelSize");
						double corr = objj.getDouble("corr");
						double med = objj.getDouble("med");
						magDbItem.addModel(modelName, trainSize, modelSize, corr, med);
					}
				}

				if (obji.has("pepLibs")) {
					JSONArray pepLibArray = obji.getJSONArray("pepLibs");
					for (int j = 0; j < pepLibArray.length(); j++) {
						JSONObject objj = pepLibArray.getJSONObject(j);
						String pepLibName = objj.getString("pepLibName");
						int pepCount = objj.getInt("pepCount");
						String protease = objj.getString("protease");
						int miss = objj.getInt("miss");
						magDbItem.addPepLib(pepLibName, pepCount, protease, miss);
					}
				}

				localList.add(magDbItem);
			}
		}

		MagDbItem[][] magDbItems = new MagDbItem[2][];
		magDbItems[0] = mgnifyList.toArray(new MagDbItem[mgnifyList.size()]);
		magDbItems[1] = localList.toArray(new MagDbItem[localList.size()]);

		return magDbItems;
	}

	public static void exportBlank(String out) {
		exportBlank(new File(out));
	}

	public static void exportBlank(File out) {
		try (PrintWriter writer = new PrintWriter(out)) {
			JSONWriter jw = new JSONWriter(writer);
			jw.object();

			jw.key(MGnify).array();

			jw.endArray();

			jw.key(localMAGs).array();
			jw.endArray();

			jw.endObject();

			writer.close();
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			LOGGER.error("Error in exporting MetaLab parameter to " + out, e);
		}
	}

	public static void export(MagDbItem[][] magDbItems, String out) {
		export(magDbItems, new File(out));
	}

	public static void export(MagDbItem[][] magDbItems, File out) {

		try (PrintWriter writer = new PrintWriter(out)) {
			JSONWriter jw = new JSONWriter(writer);
			jw.object();

			jw.key(MGnify).array();
			for (int i = 0; i < magDbItems[0].length; i++) {
				jw.object();
				jw.key("catalogueID").value(magDbItems[0][i].getCatalogueID());
				jw.key("identifier").value(magDbItems[0][i].getIdentifier());
				jw.key("repository").value(magDbItems[0][i].getRepositoryFile());

				jw.key("versionlist").array();
				ArrayList<MagVersion> versionList = magDbItems[0][i].getVersionList();
				for (int j = 0; j < versionList.size(); j++) {
					MagVersion ver = versionList.get(j);
					jw.object();
					jw.key("version").value(ver.getCatalogueVersion());
					jw.key("speciesCount").value(ver.getSpeciesCount());
					jw.key("date").value(sdf.format(ver.getDate()));
					jw.key("localAvailable").value(ver.isLocalAvailable());
					jw.key("dbName").value(ver.getDbName());
					jw.key("funcName").value(ver.getFuncName());
					jw.key("taxaName").value(ver.getTaxName());
					jw.endObject();
				}
				jw.endArray();

				jw.key("MGYG").array();
				HashMap<String, ArrayList<String>> mgygMap = magDbItems[0][i].getMgygNameMap();
				for (String parent : mgygMap.keySet()) {
					ArrayList<String> childList = mgygMap.get(parent);
					jw.object();
					jw.key("parent").value(parent);

					jw.key("children").array();
					for (String child : childList) {
						jw.object();
						jw.key("child").value(child);
						jw.endObject();
					}
					jw.endArray();
					jw.endObject();
				}
				jw.endArray();

				MagModel[] magModels = magDbItems[0][i].getMagModels();
				jw.key("models").array();
				for (MagModel magModel : magModels) {
					jw.object();
					jw.key("modelName").value(magModel.getModelName());
					jw.key("trainSize").value(magModel.getTrainSize());
					jw.key("modelSize").value(magModel.getModelSize());
					jw.key("corr").value(magModel.getCorr());
					jw.key("med").value(magModel.getMed());
					jw.endObject();
				}
				jw.endArray();

				MagPepLib[] magPepLibs = magDbItems[0][i].getMagPepLibs();
				jw.key("pepLibs").array();
				for (MagPepLib magPepLib : magPepLibs) {
					jw.object();
					jw.key("pepLibName").value(magPepLib.getLibName());
					jw.key("pepCount").value(magPepLib.getPepCount());
					jw.key("protease").value(magPepLib.getProtease());
					jw.key("miss").value(magPepLib.getMiss());
					jw.endObject();
				}
				jw.endArray();

				jw.endObject();
			}
			jw.endArray();

			jw.key(localMAGs).array();
			for (int i = 0; i < magDbItems[1].length; i++) {
				jw.object();
				jw.key("catalogueID").value(magDbItems[1][i].getCatalogueID());
				jw.key("identifier").value(magDbItems[1][i].getIdentifier());
				jw.key("repository").value(magDbItems[1][i].getRepositoryFile());

				jw.key("versionlist").array();
				ArrayList<MagVersion> versionList = magDbItems[1][i].getVersionList();
				for (int j = 0; j < versionList.size(); j++) {
					MagVersion ver = versionList.get(j);
					jw.object();
					jw.key("version").value(ver.getCatalogueVersion());
					jw.key("speciesCount").value(ver.getSpeciesCount());
					jw.key("date").value(sdf.format(ver.getDate()));
					jw.key("localAvailable").value(ver.isLocalAvailable());
					jw.key("dbName").value(ver.getDbName());
					jw.key("funcName").value(ver.getFuncName());
					jw.key("taxaName").value(ver.getTaxName());
					jw.endObject();
				}
				jw.endArray();

				MagModel[] magModels = magDbItems[1][i].getMagModels();
				jw.key("models").array();
				for (MagModel magModel : magModels) {
					jw.object();
					jw.key("modelName").value(magModel.getModelName());
					jw.key("trainSize").value(magModel.getTrainSize());
					jw.key("modelSize").value(magModel.getModelSize());
					jw.key("corr").value(magModel.getCorr());
					jw.key("med").value(magModel.getMed());
					jw.endObject();
				}
				jw.endArray();

				MagPepLib[] magPepLibs = magDbItems[1][i].getMagPepLibs();
				jw.key("pepLibs").array();
				for (MagPepLib magPepLib : magPepLibs) {
					jw.object();
					jw.key("pepLibName").value(magPepLib.getLibName());
					jw.key("pepCount").value(magPepLib.getPepCount());
					jw.key("protease").value(magPepLib.getProtease());
					jw.key("miss").value(magPepLib.getMiss());
					jw.endObject();
				}
				jw.endArray();

				jw.endObject();
			}
			jw.endArray();

			jw.endObject();

			writer.close();
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			LOGGER.error("Error in exporting MetaLab parameter to " + out, e);
		}
	}

	/**
	 * Update a single MagDbItem in the JSON configuration file This method is more
	 * efficient than exporting the entire array when only one item needs to be
	 * updated
	 * 
	 * @param magDbItem  The MagDbItem to update
	 * @param isMgnify   True if the item is from MGnify, false if it's a local MAG
	 * @param configFile The configuration file to update
	 */
	public static void updateMagDbItem(MagDbItem magDbItem, boolean isMgnify, File configFile) {
		if (!configFile.exists()) {
			LOGGER.error("Configuration file does not exist: " + configFile);
			return;
		}

		try {
			// Parse the existing configuration
			MagDbItem[][] existingItems = parse(configFile);

			// Find and update the specific item
			boolean found = false;
			if (isMgnify && existingItems[0] != null) {
				for (int i = 0; i < existingItems[0].length; i++) {
					if (existingItems[0][i].getCatalogueID().equals(magDbItem.getCatalogueID())
							&& existingItems[0][i].getIdentifier().equals(magDbItem.getIdentifier())) {
						// Update the item with new data
						existingItems[0][i] = magDbItem;
						found = true;
						break;
					}
				}
			} else if (!isMgnify && existingItems[1] != null) {
				for (int i = 0; i < existingItems[1].length; i++) {
					if (existingItems[1][i].getCatalogueID().equals(magDbItem.getCatalogueID())
							&& existingItems[1][i].getIdentifier().equals(magDbItem.getIdentifier())) {
						// Update the item with new data
						existingItems[1][i] = magDbItem;
						found = true;
						break;
					}
				}
			}

			if (!found) {
				LOGGER.warn("MagDbItem not found for update: " + magDbItem.getCatalogueID() + " - "
						+ magDbItem.getIdentifier());
				return;
			}

			// Export the updated configuration
			export(existingItems, configFile);
			LOGGER.info("Successfully updated MagDbItem: " + magDbItem.getCatalogueID() + " - "
					+ magDbItem.getIdentifier());

		} catch (Exception e) {
			LOGGER.error("Error updating MagDbItem in configuration file: " + configFile, e);
		}
	}

	/**
	 * Update a single MagDbItem in the JSON configuration file (String path
	 * version)
	 * 
	 * @param magDbItem      The MagDbItem to update
	 * @param isMgnify       True if the item is from MGnify, false if it's a local
	 *                       MAG
	 * @param configFilePath The path to the configuration file to update
	 */
	public static void updateMagDbItem(MagDbItem magDbItem, boolean isMgnify, String configFilePath) {
		updateMagDbItem(magDbItem, isMgnify, new File(configFilePath));
	}

	/**
	 * Add a new MagDbItem to the JSON configuration file
	 * 
	 * @param magDbItem  The new MagDbItem to add
	 * @param isMgnify   True if the item is from MGnify, false if it's a local MAG
	 * @param configFile The configuration file to update
	 */
	public static void addMagDbItem(MagDbItem magDbItem, boolean isMgnify, File configFile) {
		if (!configFile.exists()) {
			LOGGER.error("Configuration file does not exist: " + configFile);
			return;
		}

		try {
			// Parse the existing configuration
			MagDbItem[][] existingItems = parse(configFile);

			// Add the new item to the appropriate array
			if (isMgnify) {
				if (existingItems[0] == null) {
					existingItems[0] = new MagDbItem[] { magDbItem };
				} else {
					MagDbItem[] newArray = new MagDbItem[existingItems[0].length + 1];
					System.arraycopy(existingItems[0], 0, newArray, 0, existingItems[0].length);
					newArray[existingItems[0].length] = magDbItem;
					existingItems[0] = newArray;
				}
			} else {
				if (existingItems[1] == null) {
					existingItems[1] = new MagDbItem[] { magDbItem };
				} else {
					MagDbItem[] newArray = new MagDbItem[existingItems[1].length + 1];
					System.arraycopy(existingItems[1], 0, newArray, 0, existingItems[1].length);
					newArray[existingItems[1].length] = magDbItem;
					existingItems[1] = newArray;
				}
			}

			// Export the updated configuration
			export(existingItems, configFile);
			LOGGER.info("Successfully added new MagDbItem: " + magDbItem.getCatalogueID() + " - "
					+ magDbItem.getIdentifier());

		} catch (Exception e) {
			LOGGER.error("Error adding MagDbItem to configuration file: " + configFile, e);
		}
	}

	/**
	 * Add a new MagDbItem to the JSON configuration file (String path version)
	 * 
	 * @param magDbItem      The new MagDbItem to add
	 * @param isMgnify       True if the item is from MGnify, false if it's a local
	 *                       MAG
	 * @param configFilePath The path to the configuration file to update
	 */
	public static void addMagDbItem(MagDbItem magDbItem, boolean isMgnify, String configFilePath) {
		addMagDbItem(magDbItem, isMgnify, new File(configFilePath));
	}

}
