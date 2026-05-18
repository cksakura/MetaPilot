package uk.ac.quadram.metapilot.task.dia.par;

import uk.ac.quadram.metapilot.task.mag.MagDbItem;
import uk.ac.quadram.metapilot.task.dia.DiaRefinementMode;
import uk.ac.quadram.metapilot.task.par.MetaBaseParameter;
import uk.ac.quadram.metapilot.task.par.MetaConstants;
import uk.ac.quadram.metapilot.task.par.MetaData;

public class MetaParameterDia extends MetaBaseParameter {

	// DIA-specific fields
	private boolean librarySearch;
	private boolean selfModelSearch;
	private boolean useModel;
	private String[] library;
	private String[] models;

	private MagDbItem magDbItem;

	// Modification settings (DIA-NN / library generation)
	private boolean mExcision = true;
	private boolean carbamiC = true;
	private boolean oxidationM = false;
	private boolean acetyProN = false;
	private boolean phosphoSTY = false;
	private boolean kGG = false;
	
	// Seed workflow mode persisted for compatibility: DIRECT, FOCUS, COMPREHENSIVE, DEEP
	private String seedMode = "DIRECT";

	// Deep-seed peptide-based round mode
	private boolean usePeptideBasedRound = false;
	private double peptideDetectabilityThreshold = 0.5;
	
	public MetaParameterDia() {
		super(uk.ac.quadram.metapilot.workflow.WorkflowType.DIA_SEED.getAnalysisType(), new MetaData(), "", "",
				false, "", MetaConstants.HCD_FTMS, 1, 1, true);
		setThreadCount(Runtime.getRuntime().availableProcessors());
	}

	public MetaParameterDia(MetaBaseParameter metaPar, MagDbItem magDbItem, boolean librarySearch,
			boolean selfModelSearch, String[] library, boolean useModel, String[] models) {
		super();
		// Copy all fields from the base parameter using setters
		this.setWorkflowType(uk.ac.quadram.metapilot.workflow.WorkflowType.DIA_SEED.getAnalysisType());
		this.setMetadata(metaPar.getMetadata());
		this.setResult(metaPar.getResult());
		this.setMicroDb(metaPar.getMicroDb());
		this.setAppendHostDb(metaPar.isAppendHostDb());
		this.setHostDB(metaPar.getHostDB());
		this.setCurrentDb(metaPar.getCurrentDb());
		this.setDbSearchResultFile(metaPar.getDbSearchResultFile());
		this.setReportDir(metaPar.getReportDir());
		this.setMetafile(metaPar.getMetafile());
		this.setMs2ScanMode(metaPar.getMs2ScanMode());
		this.setCoreCount(metaPar.getCoreCount());
		this.setThreadCount(metaPar.getThreadCount());
		this.setMetaWorkflow(metaPar.isMetaWorkflow());
		this.setMscp(metaPar.getMscp());
		this.setVariMods(metaPar.getVariMods());
		this.setFixMods(metaPar.getFixMods());
		this.setEnzyme(metaPar.getEnzyme());
		this.setDigestMode(metaPar.getDigestModeString());
		this.setMissCleavages(metaPar.getMissCleavages());
		this.setReduceDb(metaPar.isReduceDb());
		this.setSpectraCluster(metaPar.isSpectraCluster());
		this.setQuantitativeMode(metaPar.getQuantitativeMode());
		this.setIsobaricTag(metaPar.getIsobaricTag());
		this.setIsobaricRefTable(metaPar.getIsobaricRefTable());
		this.setIsotopicLabels(metaPar.getIsotopicLabels());
		
		// Set DIA-specific fields
		this.magDbItem = magDbItem;
		this.librarySearch = librarySearch;
		this.selfModelSearch = selfModelSearch;
		this.useModel = useModel;
		this.library = library;
		this.models = models;
	}

	public boolean isSelfModelSearch() {
		return selfModelSearch;
	}

	public void setSelfModelSearch(boolean selfModelSearch) {
		this.selfModelSearch = selfModelSearch;
	}

	public boolean isLibrarySearch() {
		return librarySearch;
	}

	public void setLibrarySearch(boolean librarySearch) {
		this.librarySearch = librarySearch;
	}

	public boolean isUseModel() {
		return useModel;
	}

	public void setUseModel(boolean useModel) {
		this.useModel = useModel;
	}

	public String[] getLibrary() {
		return library;
	}

	public void setLibrary(String[] library) {
		this.library = library;
	}

	public String[] getModels() {
		return models;
	}

	public void setModels(String[] models) {
		this.models = models;
	}

    // Add getter and setter for MagDbItem
    public MagDbItem getMagDbItem() {
        return magDbItem;
    }
    
    public void setMagDbItem(MagDbItem magDbItem) {
        this.magDbItem = magDbItem;
        this.setMicroDb(magDbItem.getCurrentFile().getAbsolutePath());
    }

	public boolean isMExcision() { return mExcision; }
	public void setMExcision(boolean mExcision) { this.mExcision = mExcision; }
	public boolean isCarbamiC() { return carbamiC; }
	public void setCarbamiC(boolean carbamiC) { this.carbamiC = carbamiC; }
	public boolean isOxidationM() { return oxidationM; }
	public void setOxidationM(boolean oxidationM) { this.oxidationM = oxidationM; }
	public boolean isAcetyProN() { return acetyProN; }
	public void setAcetyProN(boolean acetyProN) { this.acetyProN = acetyProN; }
	public boolean isPhosphoSTY() { return phosphoSTY; }
	public void setPhosphoSTY(boolean phosphoSTY) { this.phosphoSTY = phosphoSTY; }
	public boolean isKGG() { return kGG; }
	public void setKGG(boolean kGG) { this.kGG = kGG; }
	
	public String getSeedMode() {
		return seedMode;
	}
	
	public void setSeedMode(String seedMode) {
		this.seedMode = DiaRefinementMode.fromSeedMode(seedMode).getPersistedValue();
	}

	public DiaRefinementMode getSeedRefinementMode() {
		return DiaRefinementMode.fromSeedMode(seedMode);
	}

	public void setSeedRefinementMode(DiaRefinementMode mode) {
		this.seedMode = (mode != null ? mode : DiaRefinementMode.DIRECT).getPersistedValue();
	}

	public boolean isUsePeptideBasedRound() {
		return usePeptideBasedRound;
	}

	public void setUsePeptideBasedRound(boolean usePeptideBasedRound) {
		this.usePeptideBasedRound = usePeptideBasedRound;
	}

	public double getPeptideDetectabilityThreshold() {
		return peptideDetectabilityThreshold;
	}

	public void setPeptideDetectabilityThreshold(double peptideDetectabilityThreshold) {
		this.peptideDetectabilityThreshold = peptideDetectabilityThreshold;
	}
}
