/**
 * 
 */
package uk.ac.quadram.metapilot.task.par;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintWriter;

import uk.ac.quadram.metapilot.core.mod.IsobaricTag;
import uk.ac.quadram.metapilot.task.MetaLabWorkflowType;

/**
 * 
 * MetaPilot running parameter, all the generally required parameters are
 * included in this parameter, so parameters for each version could extends this
 * class
 * 
 * @author Kai Cheng
 *
 */
public class MetaBaseParameter {

	private String workflowType;
	private MetaData metadata;
	private String result;

	private String microDb;
	private boolean appendHostDb;
	private String hostDB;
	private String currentDb;

	private File dbSearchResultFile;
	private File reportDir;
	private File metafile;

	private String ms2ScanMode;
	private int coreCount;
	private int threadCount;
	private boolean isMetaWorkflow;

	private MetaSsdbCreatePar mscp;
	private String[] variMods;
	private String[] fixMods;
	private String enzyme;
	private String digestMode;
	private int missCleavages;
	private boolean reduceDb;
	private boolean spectraCluster;

	// Quantitative mode
    private String quantitativeMode = MetaConstants.labelFree;
	private IsobaricTag isobaricTag;
	private double[][] isobaricRefTable;
	private String[][] isotopicLabels;

	private PrintWriter logWriter;

	public static final int quanResultCombine = 0;
	public static final int quanResultSeparate = 1;
	public static final int quanResultAll = 2;

	public MetaBaseParameter(String workflowType, MetaData metadata, String result, String microDb,
			boolean appendHostDb, String hostDB, String ms2ScanMode, int coreCount, int threadCount,
			boolean isMetaWorkflow) {
		this.workflowType = workflowType;
		this.metadata = metadata;
		this.result = result;
		this.microDb = microDb;
		this.appendHostDb = appendHostDb;
		this.hostDB = hostDB;
		this.ms2ScanMode = ms2ScanMode;
		this.coreCount = coreCount;
		this.threadCount = threadCount;
		this.isMetaWorkflow = isMetaWorkflow;
		this.mscp = MetaSsdbCreatePar.getDefault();
	}

	public MetaBaseParameter(String workflowType, MetaData metadata, String result, String microDb,
			boolean appendHostDb, String hostDB, String ms2ScanMode, int coreCount, int threadCount,
			boolean isMetaWorkflow, MetaSsdbCreatePar mscp) {
		this.workflowType = workflowType;
		this.metadata = metadata;
		this.result = result;
		this.microDb = microDb;
		this.appendHostDb = appendHostDb;
		this.hostDB = hostDB;
		this.ms2ScanMode = ms2ScanMode;
		this.coreCount = coreCount;
		this.threadCount = threadCount;
		this.isMetaWorkflow = isMetaWorkflow;
	}

	public String[][] getIsotopicLabels() {
		return isotopicLabels;
	}

	public void setIsotopicLabels(String[][] isotopicLabels) {
		this.isotopicLabels = isotopicLabels;
	}

	public int getMissCleavages() {
		return missCleavages;
	}

	public void setMissCleavages(int missCleavages) {
		this.missCleavages = missCleavages;
	}

	public boolean isReduceDb() {
		return reduceDb;
	}

	public void setReduceDb(boolean reduceDb) {
		this.reduceDb = reduceDb;
	}

	public boolean isSpectraCluster() {
		return spectraCluster;
	}

	public void setSpectraCluster(boolean spectraCluster) {
		this.spectraCluster = spectraCluster;
	}

	public String[] getVariMods() {
		return variMods;
	}

	public void setVariMods(String[] variMods) {
		this.variMods = variMods;
	}

	public String[] getFixMods() {
		return fixMods;
	}

	public void setFixMods(String[] fixMods) {
		this.fixMods = fixMods;
	}

	public String getEnzyme() {
		return enzyme;
	}

	public void setEnzyme(String enzyme) {
		this.enzyme = enzyme;
	}
	
	public String getDigestModeString() {
		return digestMode;
	}

	public int getDigestMode() {
		if (digestMode.equals("Semi specific")) {
			return 1;
		} else if (digestMode.equals("Unspecific")) {
			return 0;
		} else {
			return 3;
		}
	}

	public void setDigestMode(String digestMode) {
		this.digestMode = digestMode;
	}

	public String getQuantitativeMode() {
		return quantitativeMode;
	}

	public void setQuantitativeMode(String quantitativeMode) {
		this.quantitativeMode = quantitativeMode;
	}

	public IsobaricTag getIsobaricTag() {
		return isobaricTag;
	}

	public void setIsobaricTag(IsobaricTag isobaricTag) {
		this.isobaricTag = isobaricTag;
	}

	public double[][] getIsobaricRefTable() {
		return isobaricRefTable;
	}

	public void setIsobaricRefTable(double[][] isobaricRefTable) {
		this.isobaricRefTable = isobaricRefTable;
	}

	public void setLogWriter(PrintWriter logWriter) {
		this.logWriter = logWriter;
	}

	public MetaData getMetadata() {
		return metadata;
	}

	public void setMetadata(MetaData metadata) {
		this.metadata = metadata;
	}

	public File getMetafile() {
		return metafile;
	}

	public void setMetafile(File metafile) {
		this.metafile = metafile;
	}

	public int getCoreCount() {
		return coreCount;
	}

	public void setCoreCount(int coreCount) {
		this.coreCount = coreCount;
	}

	public int getThreadCount() {
		return threadCount;
	}

	public void setThreadCount(int threadCount) {
		this.threadCount = threadCount;
	}

	public String getMs2ScanMode() {
		return ms2ScanMode;
	}

	public void setMs2ScanMode(String ms2ScanMode) {
		this.ms2ScanMode = ms2ScanMode;
	}

	/**
	 * the result path is a folder containing all of the result files
	 * 
	 * @return the directory of the result file
	 */
	public String getResult() {
		return result;
	}

	public void setResult(String result) {
		this.result = result;
	}

	public String getMicroDb() {
		return microDb;
	}

	public void setMicroDb(String microDb) {
		this.microDb = microDb;
	}

	public boolean isAppendHostDb() {
		return appendHostDb;
	}

	public void setAppendHostDb(boolean appendHostDb) {
		this.appendHostDb = appendHostDb;
	}

	public String getHostDB() {
		return hostDB;
	}

	public void setHostDB(String hostDB) {
		this.hostDB = hostDB;
	}

	public String getCurrentDb() {
		if (currentDb == null) {
			currentDb = new String(this.microDb);
		}
		return currentDb;
	}

	public void setCurrentDb(String currentDb) {
		this.currentDb = currentDb;
	}

	public MetaLabWorkflowType getWorkflowType() {
		return MetaLabWorkflowType.valueOf(workflowType);
	}

	public void setWorkflowType(String workflowType) {
		this.workflowType = workflowType;
	}

	/**
	 * If true, do the taxonomy analysis and functional annotation
	 * 
	 * @return
	 */
	public boolean isMetaWorkflow() {
		return isMetaWorkflow;
	}

	public void setMetaWorkflow(boolean isMetaWorkflow) {
		this.isMetaWorkflow = isMetaWorkflow;
	}

	public File getDbSearchResultFile() {
		return dbSearchResultFile;
	}

	public void setDbSearchResultFile(File dbSearchResultFile) {
		this.dbSearchResultFile = dbSearchResultFile;
	}

	public File getReportDir() {
		return reportDir;
	}

	public void setReportDir(File reportDir) {
		this.reportDir = reportDir;
	}

	public File getSpectraFile() {
		File spectraFile = new File(result, "spectra");
		if (!spectraFile.exists()) {
			spectraFile.mkdir();
		}
		return spectraFile;
	}

	public MetaSsdbCreatePar getMscp() {
		return mscp;
	}

	public void setMscp(MetaSsdbCreatePar mscp) {
		this.mscp = mscp;
	}

	public PrintWriter getLogWriter() throws FileNotFoundException {
		if (this.logWriter == null) {
			logWriter = new PrintWriter(new File(result, "log.txt"));
		}
		return this.logWriter;
	}

}
