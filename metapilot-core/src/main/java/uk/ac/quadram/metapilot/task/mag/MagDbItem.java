package uk.ac.quadram.metapilot.task.mag;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.Files;
import java.nio.file.StandardCopyOption;
import java.util.ArrayList;
import java.util.Date;
import java.util.HashMap;
import java.util.HashSet;

import org.apache.commons.io.FileUtils;

import uk.ac.quadram.metapilot.dbSearch.diann.DiaNNParquetReader;

public class MagDbItem {

	public static final String MARKER2 = "Marker2";
	public static final String MARKER4 = "Marker4";
	/** Subdirectory of a MAG database folder holding named spectral libraries (DIA Seed, etc.). */
	public static final String LIBRARIES_FOLDER_NAME = "libraries";

	private String catalogueID;
	private String identifier;
	private MagVersion version;
	private ArrayList<MagVersion> versionList;
	private HashMap<String, ArrayList<String>> mgygNameMap;
	private MagVersion usedVersion;

	private File repositoryFile;
	private File currentFile;
	private File hostDbFile;
	
	private HashMap<String, MagModel> magModelMap;
	private HashMap<String, MagPepLib> magPepLibMap;

	public MagDbItem(String catalogueID, String identifier) {
		this.catalogueID = catalogueID;
		this.identifier = identifier;
		this.versionList = new ArrayList<MagVersion>();

		this.repositoryFile = null;
		
		this.magModelMap = new HashMap<String, MagModel>();
		this.magPepLibMap = new HashMap<String, MagPepLib>();
	}

	public MagDbItem(String catalogueID, String identifier, File directoryFile) {
		this.catalogueID = catalogueID;
		this.identifier = identifier;
		this.versionList = new ArrayList<MagVersion>();

		this.repositoryFile = directoryFile;
		
		this.magModelMap = new HashMap<String, MagModel>();
		this.magPepLibMap = new HashMap<String, MagPepLib>();
	}

	public void addVersion(String catalogueVersion, int speciesCount, Date date, boolean localAvailable) {
		MagVersion version = new MagVersion(catalogueVersion, speciesCount, date, localAvailable);
		this.versionList.add(version);
		this.version = version;
	}
	
	public void addVersion(String catalogueVersion, int speciesCount, Date date, boolean localAvailable, String dbName,
			String funcName, String taxaName) {
		MagVersion version = new MagVersion(catalogueVersion, speciesCount, date, localAvailable, dbName, funcName,
				taxaName);
		this.versionList.add(version);
		this.version = version;
	}

	public MagModel[] getMagModels() {
		MagModel[] magModels = this.magModelMap.values().toArray(MagModel[]::new);
		return magModels;
	}

	public MagPepLib[] getMagPepLibs() {
		MagPepLib[] magPepLibs = this.magPepLibMap.values().toArray(MagPepLib[]::new);
		return magPepLibs;
	}

	public void addModel(String modelName, int trainSize, int modelSize, double corr, double med) {
		this.magModelMap.put(modelName, new MagModel(modelName, trainSize, modelSize, corr, med));
	}

	public void addPepLib(String pepLibName, int pepCount, String protease, int miss) {
		this.magPepLibMap.put(pepLibName, new MagPepLib(pepLibName, pepCount, protease, miss, true));
	}

	public void addPepLib(String pepLibName, int pepCount, String protease, int miss, boolean matchToMag) {
		this.magPepLibMap.put(pepLibName, new MagPepLib(pepLibName, pepCount, protease, miss, matchToMag));
	}
	
	public boolean removeModel(MagModel magModel) {
		String modelName = magModel.getModelName();
		MagModel removedModel = this.magModelMap.remove(modelName);
		return removedModel != null;
	}
	
	public boolean removePepLib(MagPepLib magPepLib) {
		String libName = magPepLib.getLibName();
		MagPepLib removedPepLib = this.magPepLibMap.remove(libName);
		return removedPepLib != null;
	}
	
	public void removeAllModels() {
		this.magModelMap.clear();
	}
	
	public void removeAllPepLibs() {
		this.magPepLibMap.clear();
	}
	
	public void addModelFromFile(File file) {
		String modelName = file.getName();
		File infoFile = new File(file, modelName + "_info.txt");
		if (infoFile.exists()) {
			try (BufferedReader reader = new BufferedReader(new FileReader(infoFile))) {
				String line = reader.readLine();
				String[] cs = line.split("\t");
				int trainSize = Integer.parseInt(cs[1]);

				line = reader.readLine();
				cs = line.split("\t");
				int dataSize = Integer.parseInt(cs[1]);

				line = reader.readLine();
				cs = line.split("\t");
				double corr = Double.parseDouble(cs[1]);

				line = reader.readLine();
				cs = line.split("\t");
				double med = Double.parseDouble(cs[1]);

				addModel(modelName, trainSize, dataSize, corr, med);

				reader.close();
			} catch (Exception e) {
				System.err.println(e);
			}
		}
	}
	
	public void addPepLibFromFile(File file) {
		String pepLibName = file.getName();
		File infoFile = new File(file, pepLibName + "_info.txt");
		if (infoFile.exists()) {
			try (BufferedReader reader = new BufferedReader(new FileReader(infoFile))) {
				String line;
				int pepCount = 0;
				String enzyme = "Trypsin/P";
				int miss = 1;
				boolean matchToMag = true;
				while ((line = reader.readLine()) != null) {
					String[] cs = line.split("\t");
					if (cs.length < 2) continue;
					String key = cs[0].trim().replaceAll(":$", "");
					if (key.equalsIgnoreCase("Data set size")) pepCount = Integer.parseInt(cs[1].trim());
					else if (key.equalsIgnoreCase("Enzyme")) enzyme = cs[1].trim();
					else if (key.equalsIgnoreCase("Miss cleavages")) miss = Integer.parseInt(cs[1].trim());
					else if (key.equalsIgnoreCase("MatchToMag")) matchToMag = Boolean.parseBoolean(cs[1].trim());
				}
				addPepLib(pepLibName, pepCount, enzyme, miss, matchToMag);
			} catch (Exception e) {
				System.err.println(e);
			}
		}
	}
	
	public File getPredictedPepFile() {
		if (currentFile == null) {
			return null;
		}

		File predictedPepFile = new File(currentFile, "predicted");
		return predictedPepFile;
	}
	
	public void resetMgygMap() {
		mgygNameMap = new HashMap<String, ArrayList<String>>();
	}
	
	public void addMgygName(String parent, String child) {
		if (mgygNameMap == null) {
			mgygNameMap = new HashMap<String, ArrayList<String>>();
		}

		if (!mgygNameMap.containsKey(parent)) {
			mgygNameMap.put(parent, new ArrayList<String>());
		}

		mgygNameMap.get(parent).add(child);
	}

	public HashMap<String, ArrayList<String>> getMgygNameMap() {
		return mgygNameMap;
	}

	public String getCatalogueID() {
		return catalogueID;
	}

	public void setCatalogueID(String catalogueID) {
		this.catalogueID = catalogueID;
	}

	public ArrayList<MagVersion> getVersionList() {
		return versionList;
	}

	public String getCatalogueVersion() {
		if (version == null) {
			return "";
		}
		return version.getCatalogueVersion();
	}

	public int getSpeciesCount() {
		if (version == null) {
			return 0;
		}
		return version.getSpeciesCount();
	}
	
	public void setAvailable() {
		if (version != null) {
			version.localAvailable = true;
		}
	}
	
	public boolean isAvailable() {
		if (version == null) {
			return false;
		}
		return version.isLocalAvailable();
	}
	
	public String getOldVersion() {
		if (versionList.size() == 0) {
			return "";
		} else if (versionList.size() == 1) {
			return version.getCatalogueVersion();
		} else {
			return versionList.get(versionList.size() - 2).getCatalogueVersion();
		}
	}

	public String getIdentifier() {
		return identifier;
	}

	public void setIdentifier(String identifier) {
		this.identifier = identifier;
	}

	public Date getDate() {
		if (version == null) {
			return null;
		}
		return version.getDate();
	}

	public MagVersion getUsedVersion() {
		if (usedVersion == null) {
			for (int i = 0; i < versionList.size(); i++) {
				MagVersion veri = versionList.get(i);
				if (veri.isLocalAvailable()) {
					usedVersion = veri;
				}
			}
		}
		return usedVersion;
	}

	public void setUsedVersion(MagVersion usedVersion) {
		this.usedVersion = usedVersion;
		this.currentFile = new File(repositoryFile, usedVersion.getCatalogueVersion());
		if (!this.currentFile.exists()) {
			this.currentFile = repositoryFile;
		}
	}

	public File getRepositoryFile() {
		return repositoryFile;
	}

	public void setRepositoryFile(File repositoryFile) {
		this.repositoryFile = repositoryFile;
	}

	public boolean updateAvailable() {
		if (repositoryFile == null || !repositoryFile.exists()) {
			return true;
		}

		if (version == null) {
			return true;
		}

		File currentFolder = new File(repositoryFile, this.version.getCatalogueVersion());
		if (!currentFolder.exists()) {
			return true;
		}
		return false;
	}

	public File getHostDbFile() {
		return hostDbFile;
	}

	public void setHostDbFile(File hostDbFile) {
		this.hostDbFile = hostDbFile;
	}

	public File getCurrentFile() {
		return currentFile;
	}

	public void setCurrentFile(File currentFile) {
		this.currentFile = currentFile;
	}

	public File getHapFastaFile() {
		return getRibElonFastaFile();
	}

	public File getRibElonFastaFile() {
		if (repositoryFile == null || !repositoryFile.exists()) {
			return null;
		}
		if (currentFile == null || !currentFile.exists()) {
			return null;
		}

		File marker2File = new File(currentFile, "Marker2.fasta");
		if (!marker2File.exists()) {
			File ribFile = new File(currentFile, "rib_elon.fasta");
			if (ribFile.exists()) {
				ribFile.renameTo(marker2File);
			}
		}

		return marker2File;
	}

	/**
	 * Returns the path to Marker4_highpep.fasta (used when caller knows it exists).
	 */
	public File getMarkerHighPepFastaFile() {
		if (repositoryFile == null || !repositoryFile.exists()) {
			return null;
		}
		if (currentFile == null || !currentFile.exists()) {
			return null;
		}
		return new File(currentFile, "Marker4_highpep.fasta");
	}

	/**
	 * Returns the effective Marker4 FASTA for library generation: Marker4_highpep.fasta if it exists
	 * (generated when Marker4.fasta has &gt;100k proteins), otherwise Marker4.fasta.
	 */
	public File getMarker4FastaForLibrary() {
		if (repositoryFile == null || !repositoryFile.exists()) {
			return null;
		}
		if (currentFile == null || !currentFile.exists()) {
			return null;
		}
		File highpep = new File(currentFile, "Marker4_highpep.fasta");
		if (highpep.exists()) {
			return highpep;
		}
		return new File(currentFile, "Marker4.fasta");
	}

	/**
	 * {@code libraries} under the current MAG database folder, or null if unavailable.
	 */
	public File getLibrariesRootFolder() {
		if (repositoryFile == null || !repositoryFile.exists()) {
			return null;
		}
		if (currentFile == null || !currentFile.exists()) {
			return null;
		}
		return new File(currentFile, LIBRARIES_FOLDER_NAME);
	}

	/**
	 * Canonical folder for the Marker4 spectral library ({@code libraries/Marker4}), discoverable as library id
	 * {@link #MARKER4} by DIA Seed.
	 */
	public File getMarker4LibraryFolder() {
		File libRoot = getLibrariesRootFolder();
		if (libRoot == null) {
			return null;
		}
		return new File(libRoot, MARKER4);
	}

	/**
	 * Canonical Marker4 spectral library path: {@code libraries/Marker4/Marker4.speclib}. Used when generating
	 * the library and after migration from legacy root-level files.
	 */
	public File getMarker4LibFile() {
		File folder = getMarker4LibraryFolder();
		if (folder == null) {
			return null;
		}
		return new File(folder, MARKER4 + ".speclib");
	}

	/**
	 * Creates {@code libraries/Marker4} (and {@code libraries} if needed) before writing the speclib, parquet, or
	 * info files. Legacy move also mkdirs, but fresh generation from FASTA must call this first.
	 *
	 * @return false if paths are invalid or the directory could not be created
	 */
	public boolean ensureMarker4LibraryFolderExists() {
		File folder = getMarker4LibraryFolder();
		if (folder == null) {
			return false;
		}
		try {
			FileUtils.forceMkdir(folder);
			return folder.isDirectory();
		} catch (IOException e) {
			return false;
		}
	}

	private File getMarker4LibLegacyPredictedFile() {
		if (repositoryFile == null || !repositoryFile.exists()) {
			return null;
		}
		if (currentFile == null || !currentFile.exists()) {
			return null;
		}
		File highpep = new File(currentFile, "Marker4_highpep.fasta");
		if (highpep.exists()) {
			return new File(currentFile, "Marker4_highpep.predicted.speclib");
		}
		return new File(currentFile, "Marker4.predicted.speclib");
	}

	/**
	 * Removes legacy root-level predicted speclibs when the canonical library already lives under {@code libraries/Marker4}.
	 */
	private void deleteLegacyMarker4PredictedSpeclibsInParent() {
		if (currentFile == null || !currentFile.exists()) {
			return;
		}
		FileUtils.deleteQuietly(new File(currentFile, "Marker4.predicted.speclib"));
		FileUtils.deleteQuietly(new File(currentFile, "Marker4_highpep.predicted.speclib"));
	}

	/**
	 * If {@code libraries/Marker4/Marker4.speclib} already exists, deletes any legacy root-level
	 * {@code Marker4*.predicted.speclib}. Otherwise, if a legacy predicted speclib exists in the MAG folder root,
	 * moves it into {@code libraries/Marker4/Marker4.speclib} (previous MetaPilot layout).
	 *
	 * @return true if {@link #getMarker4LibFile()} exists after the call
	 */
	public boolean moveLegacyMarker4SpeclibToCanonical() {
		File canonical = getMarker4LibFile();
		if (canonical == null) {
			return false;
		}
		if (canonical.exists()) {
			deleteLegacyMarker4PredictedSpeclibsInParent();
			return true;
		}
		File legacy = getMarker4LibLegacyPredictedFile();
		if (legacy == null || !legacy.exists()) {
			return false;
		}
		try {
			FileUtils.forceMkdir(canonical.getParentFile());
			Files.move(legacy.toPath(), canonical.toPath(), StandardCopyOption.REPLACE_EXISTING);
		} catch (IOException e) {
			return canonical.exists();
		}
		return canonical.exists();
	}

	/**
	 * Resolves an existing Marker4 spectral library for DiaNN: canonical {@code libraries/Marker4} (.speclib,
	 * .tsv, .parquet) first, then legacy root-level {@code Marker4*.predicted.speclib}.
	 *
	 * @return an existing file, or null
	 */
	public File resolveMarker4LibraryFile() {
		File canonicalFolder = getMarker4LibraryFolder();
		if (canonicalFolder != null && canonicalFolder.isDirectory()) {
			for (String ext : new String[] { ".speclib", ".tsv", ".parquet" }) {
				File f = new File(canonicalFolder, MARKER4 + ext);
				if (f.exists()) {
					return f;
				}
			}
		}
		File legacy = getMarker4LibLegacyPredictedFile();
		if (legacy != null && legacy.exists()) {
			return legacy;
		}
		return null;
	}

	/**
	 * Resolves a spectral library id under this MAG database to an existing file. Callers should treat an absolute
	 * path first ({@code new File(name).exists()}). For {@link #MARKER4}, applies canonical layout and legacy
	 * fallback ({@link #resolveMarker4LibraryFile()}).
	 */
	public File resolveSpectralLibraryFile(String libraryName) {
		if (libraryName == null || libraryName.isEmpty()) {
			return null;
		}
		if (MARKER4.equals(libraryName)) {
			moveLegacyMarker4SpeclibToCanonical();
			return resolveMarker4LibraryFile();
		}
		if (repositoryFile == null || !repositoryFile.exists()) {
			return null;
		}
		if (currentFile == null || !currentFile.exists()) {
			return null;
		}
		File libFolder = new File(new File(currentFile, LIBRARIES_FOLDER_NAME), libraryName);
		for (String ext : new String[] { ".speclib", ".tsv", ".parquet" }) {
			File f = new File(libFolder, libraryName + ext);
			if (f.exists()) {
				return f;
			}
		}
		return null;
	}

	/**
	 * Moves legacy Marker4 speclib into {@code libraries/Marker4} if needed, then refreshes {@link #getMagPepLibs()}
	 * from {@code Marker4_info.txt} or {@code Marker4.parquet} when present. The info file and parquet are produced
	 * when a Marker DIA task completes spectral library setup in the GUI module.
	 */
	public void registerMarker4SpectralLibraryIfPresent() {
		moveLegacyMarker4SpeclibToCanonical();
		File folder = getMarker4LibraryFolder();
		if (folder == null || !folder.isDirectory()) {
			return;
		}
		File infoFile = new File(folder, MARKER4 + "_info.txt");
		if (infoFile.exists()) {
			addPepLibFromFile(folder);
			return;
		}
		File parquet = new File(folder, MARKER4 + ".parquet");
		if (parquet.exists()) {
			writeMarker4InfoAndRegisterFromParquet(parquet, folder);
		}
	}

	private void writeMarker4InfoAndRegisterFromParquet(File parquet, File libraryFolder) {
		try {
			HashSet<String> pepSet = DiaNNParquetReader.getPeptideSet(parquet.getAbsolutePath());
			int n = pepSet != null ? pepSet.size() : 0;
			File infoFile = new File(libraryFolder, MARKER4 + "_info.txt");
			try (PrintWriter writer = new PrintWriter(infoFile)) {
				writer.println("Data set size:\t" + n);
				writer.println("Enzyme:\tTrypsin/P");
				writer.println("Miss cleavages:\t1");
				writer.println("MatchToMag:\ttrue");
			}
			addPepLib(MARKER4, n, "Trypsin/P", 1, true);
		} catch (Exception e) {
			System.err.println("Failed to sync Marker4 library from parquet: " + e.getMessage());
		}
	}

	public File getHapLibFile(String selectorMode) {
		return getMarker4LibFile();
	}

	/** Delegates to {@link #getMarker4LibFile()}. Kept for backward compatibility. */
	public File getMarkerHighPepLibFile() {
		return getMarker4LibFile();
	}

	public File getHapFastaFile(String selectorMode) {
		if (MARKER4.equalsIgnoreCase(selectorMode)) {
			return getMarker4FastaForLibrary();
		}
		return getRibElonFastaFile();
	}
	
	public File getHapFastaHdfFile() {
		if (repositoryFile == null || !repositoryFile.exists()) {
			return null;
		}
		if (currentFile == null || !currentFile.exists()) {
			return null;
		}

		return new File(currentFile, "Marker2.hdf");
	}

	public File getHapFastaTDFile() {
		return getHapFastaTDFile(null);
	}

	/**
	 * Returns the target-decoy FASTA path for the given genome selector mode.
	 * Marker4 uses Marker4_td.fasta; otherwise Marker2_td.fasta.
	 * The file may not exist yet (caller should create it via FastaManager.addDecoy if needed).
	 */
	public File getHapFastaTDFile(String selectorMode) {
		if (repositoryFile == null || !repositoryFile.exists()) {
			return null;
		}
		if (currentFile == null || !currentFile.exists()) {
			return null;
		}
		if (MARKER4.equalsIgnoreCase(selectorMode)) {
			return new File(currentFile, "Marker4_td.fasta");
		}
		return new File(currentFile, "Marker2_td.fasta");
	}

	public File getOriginalDbFolder() {
		if (repositoryFile == null || !repositoryFile.exists()) {
			return null;
		}
		if (usedVersion == null) {
			return null;
		}

		return new File(currentFile, usedVersion.getDbName());
	}
	
	public File getOriginalFuncFolder() {
		if (repositoryFile == null || !repositoryFile.exists()) {
			return null;
		}
		if (usedVersion == null) {
			return null;
		}

		return new File(currentFile, usedVersion.getFuncName());
	}
	
	public File getOriginalTaxaFile() {
		if (repositoryFile == null || !repositoryFile.exists()) {
			return null;
		}
		if (usedVersion == null) {
			return null;
		}

		return new File(currentFile, usedVersion.getTaxName());
	}
	
	public class MagPepLib {
		private String libName;
		private int pepCount;
		private String protease;
		private int miss;
		/** True = microbiome (match peptides to MAG for protein/genome info); false = standalone (host/annotated, no MAG matching). */
		private boolean matchToMag;

		public MagPepLib(String libName, int pepCount, String protease, int miss) {
			this(libName, pepCount, protease, miss, true);
		}

		public MagPepLib(String libName, int pepCount, String protease, int miss, boolean matchToMag) {
			this.libName = libName;
			this.pepCount = pepCount;
			this.protease = protease;
			this.miss = miss;
			this.matchToMag = matchToMag;
		}

		public String getLibName() {
			return libName;
		}

		public int getPepCount() {
			return pepCount;
		}

		public String getProtease() {
			return protease;
		}

		public int getMiss() {
			return miss;
		}

		public boolean isMatchToMag() {
			return matchToMag;
		}

		public void setMatchToMag(boolean matchToMag) {
			this.matchToMag = matchToMag;
		}

		public void delete() {
			File modelFolder = new File(new File(currentFile, "libraries"), libName);
			FileUtils.deleteQuietly(modelFolder);
		}

	}
	
	public class MagModel {
		private String modelName;
		private int trainSize;
		private int modelSize;
		private double corr;
		private double med;

		public MagModel(String modelName, int trainSize, int modelSize, double corr, double med) {
			this.modelName = modelName;
			this.trainSize = trainSize;
			this.modelSize = modelSize;
			this.corr = corr;
			this.med = med;
		}

		public String getModelName() {
			return modelName;
		}

		public void setModelName(String modelName) {
			this.modelName = modelName;
		}

		public int getTrainSize() {
			return trainSize;
		}

		public void setTrainSize(int trainSize) {
			this.trainSize = trainSize;
		}

		public int getModelSize() {
			return modelSize;
		}

		public void setModelSize(int modelSize) {
			this.modelSize = modelSize;
		}

		public double getCorr() {
			return corr;
		}

		public void setCorr(double corr) {
			this.corr = corr;
		}

		public double getMed() {
			return med;
		}

		public void setMed(double med) {
			this.med = med;
		}

		public File getModelFile() {
			File modelFolder = new File(new File(currentFile, "models"), modelName);
			if (!modelFolder.exists()) {
				return null;
			}

			File modelFile = new File(modelFolder, modelName + ".h5");
			if (!modelFile.exists()) {
				return null;
			}

			return modelFile;
		}

		public void delete() {
			File modelFolder = new File(new File(currentFile, "models"), modelName);
			FileUtils.deleteQuietly(modelFolder);
		}
	}
	
	public class MagVersion {

		private String catalogueVersion;
		private int speciesCount;
		private Date date;
		private boolean localAvailable;
		
		private String dbName;
		private String funcName;
		private String taxName;

		MagVersion(String catalogueVersion, int speciesCount, Date date, boolean localAvailable) {
			this.catalogueVersion = catalogueVersion;
			this.speciesCount = speciesCount;
			this.date = date;
			this.localAvailable = localAvailable;
			this.dbName = "original_db";
			this.funcName = "eggNOG";
			this.taxName = "genomes-all_metadata.tsv";
		}
		
		MagVersion(String catalogueVersion, int speciesCount, Date date, boolean localAvailable, String dbName, String funcName, String taxName) {
			this.catalogueVersion = catalogueVersion;
			this.speciesCount = speciesCount;
			this.date = date;
			this.localAvailable = localAvailable;
			this.dbName = dbName;
			this.funcName = funcName;
			this.taxName = taxName;
		}

		public String getCatalogueVersion() {
			return catalogueVersion;
		}

		public int getSpeciesCount() {
			return speciesCount;
		}

		public Date getDate() {
			return date;
		}

		public boolean isLocalAvailable() {
			return localAvailable;
		}
		
		public String getDbName() {
			return dbName;
		}

		public String getFuncName() {
			return funcName;
		}

		public String getTaxName() {
			return taxName;
		}

		public String toString() {
			return catalogueVersion;
		}
	}

	public String toString() {
		return catalogueID;
	}

	public int hashCode() {
		return (catalogueID + "_" + this.usedVersion.catalogueVersion).hashCode();
	}

	public boolean equals(Object obj) {
		if (obj instanceof MagDbItem) {
			MagDbItem magDbItem = (MagDbItem) obj;
			if (!magDbItem.catalogueID.equals(this.catalogueID)) {
				return false;
			}

			if (magDbItem.usedVersion == null || this.usedVersion == null) {
				return true;
			}

			if (!magDbItem.usedVersion.catalogueVersion.equals(this.usedVersion.catalogueVersion)) {
				return false;
			}
			return true;
		} else {
			return false;

		}
	}

	public Object[] getTableRow() {
		Object[] objs = new Object[] { isAvailable(), catalogueID, getCatalogueVersion(), getSpeciesCount(),
				updateAvailable() };
		return objs;
	}

	public static Object[] getColumnNames() {
		Object[] objs = new Object[] { "Saved in local", "Catalogue ID", "Version", "Species Count",
				"Update available" };
		return objs;
	}
}
