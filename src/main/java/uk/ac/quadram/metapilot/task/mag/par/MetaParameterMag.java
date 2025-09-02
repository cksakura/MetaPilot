package uk.ac.quadram.metapilot.task.mag.par;

import java.io.File;
import java.util.ArrayList;

import uk.ac.quadram.metapilot.core.mod.IsobaricTag;
import uk.ac.quadram.metapilot.dbSearch.MaxQuant.MaxquantModification;
import uk.ac.quadram.metapilot.task.mag.MagDbItem;
import uk.ac.quadram.metapilot.task.mag.MagDbItem.MagVersion;
import uk.ac.quadram.metapilot.task.par.MetaData;
import uk.ac.quadram.metapilot.task.par.MetaBaseParameter;
import uk.ac.quadram.metapilot.task.par.MetaPep2TaxaPar;
import uk.ac.quadram.metapilot.task.v2.par.MetaOpenPar;

public class MetaParameterMag extends MetaBaseParameter {

	// Fields from MetaParameterPFind
	private String[] fixMods;
	private String[] variMods;
	private String enzyme;
	private int missCleavages;
	private int digestMode;
	private String quanMode;
	private MaxquantModification[] isobaric;
	private IsobaricTag isobaricTag;
	private boolean combineLabel;
	private double[][] isoCorFactor;
	protected MetaOpenPar osp;
	protected MetaPep2TaxaPar ptp;

	// Fields specific to MetaParameterMag
	protected MagDbItem[][] magDbItems;
	protected String[] fragpipeEnzyme;
	protected MagDbItem usedMagDbItem;
	protected String magDb;
	protected String magDbVersion;

	protected MaxquantModification[] fullFixMods;
	protected MaxquantModification[] fullVariMods;

	public MetaParameterMag(MetaBaseParameter metaPar) {
		super(metaPar.getWorkflowType().name(), metaPar.getMetadata(), metaPar.getResult(), 
			  metaPar.getMicroDb(), metaPar.isAppendHostDb(), metaPar.getHostDB(), 
			  metaPar.getMs2ScanMode(), metaPar.getCoreCount(), metaPar.getThreadCount(), 
			  metaPar.isMetaWorkflow(), metaPar.getMscp());
		
		// Initialize inherited fields from MetaParameterPFind with defaults
		this.fixMods = new String[] { "Carbamidomethyl[C] 0" };
		this.variMods = new String[] { "Acetyl[ProteinN-term] 0", "Oxidation[M] 0" };
		this.enzyme = "Trypsin_P";
		this.missCleavages = 2;
		this.digestMode = 0;
		this.quanMode = "Label-free";
		this.isobaric = new MaxquantModification[] {};
		this.isobaricTag = null;
		this.combineLabel = false;
		this.isoCorFactor = new double[][] {};
		this.osp = MetaOpenPar.getDefault();
		this.ptp = MetaPep2TaxaPar.getDefault();
	}

	public MetaParameterMag(MetaBaseParameter metaPar, String[] fixMods, String[] variMods, String enzyme,
			int missCleavages, int digestMode, String quanMode, MaxquantModification[] isobaric,
			IsobaricTag isobaricTag, boolean combineLabel, double[][] isoCorFactor, String magDb, String magDbVersion) {
		super(metaPar.getWorkflowType().name(), metaPar.getMetadata(), metaPar.getResult(), 
			  metaPar.getMicroDb(), metaPar.isAppendHostDb(), metaPar.getHostDB(), 
			  metaPar.getMs2ScanMode(), metaPar.getCoreCount(), metaPar.getThreadCount(), 
			  metaPar.isMetaWorkflow(), metaPar.getMscp());
		
		// Initialize inherited fields from MetaParameterPFind
		this.fixMods = fixMods;
		this.variMods = variMods;
		this.enzyme = enzyme;
		this.missCleavages = missCleavages;
		this.digestMode = digestMode;
		this.quanMode = quanMode;
		this.isobaric = isobaric;
		this.isobaricTag = isobaricTag;
		this.combineLabel = combineLabel;
		this.isoCorFactor = isoCorFactor;
		this.osp = MetaOpenPar.getDefault();
		this.ptp = MetaPep2TaxaPar.getDefault();
		
		// Initialize MAG-specific fields
		this.magDb = magDb;
		this.magDbVersion = magDbVersion;
	}

	public MetaParameterMag(MetaBaseParameter metaPar, String[] fixMods, String[] variMods, String enzyme,
			int missCleavages, int digestMode, String quanMode, MaxquantModification[] isobaric,
			IsobaricTag isobaricTag, boolean combineLabel, double[][] isoCorFactor, String magDb, String magDbVersion,
			MagDbItem[][] magDbItems) {
		super(metaPar.getWorkflowType().name(), metaPar.getMetadata(), metaPar.getResult(), 
			  metaPar.getMicroDb(), metaPar.isAppendHostDb(), metaPar.getHostDB(), 
			  metaPar.getMs2ScanMode(), metaPar.getCoreCount(), metaPar.getThreadCount(), 
			  metaPar.isMetaWorkflow(), metaPar.getMscp());
		
		// Initialize inherited fields from MetaParameterPFind
		this.fixMods = fixMods;
		this.variMods = variMods;
		this.enzyme = enzyme;
		this.missCleavages = missCleavages;
		this.digestMode = digestMode;
		this.quanMode = quanMode;
		this.isobaric = isobaric;
		this.isobaricTag = isobaricTag;
		this.combineLabel = combineLabel;
		this.isoCorFactor = isoCorFactor;
		this.osp = MetaOpenPar.getDefault();
		this.ptp = MetaPep2TaxaPar.getDefault();
		
		// Initialize MAG-specific fields
		this.magDb = magDb;
		this.magDbVersion = magDbVersion;
		this.magDbItems = magDbItems;
		
		// Find and set the used MAG database item
		for (int i = 0; i < magDbItems.length; i++) {
			for (int j = 0; j < magDbItems[i].length; j++) {
				if (magDb.equals(magDbItems[i][j].getCatalogueID())) {
					usedMagDbItem = magDbItems[i][j];

					ArrayList<MagVersion> verlist = magDbItems[i][j].getVersionList();
					for (int k = 0; k < verlist.size(); k++) {
						if (verlist.get(k).getCatalogueVersion().equals(magDbVersion)) {
							usedMagDbItem.setUsedVersion(verlist.get(k));
						}
					}
				}
			}
		}
	}

	public boolean hasRawFiles() {
		MetaData metaData = this.getMetadata();
		if (metaData == null) {
			return false;
		} else {
			String[] raws = metaData.getRawFiles();
			if (raws == null || raws.length == 0) {
				return false;
			} else {
				for (int i = 0; i < raws.length; i++) {
					if (!(new File(raws[i])).exists()) {
						return false;
					}
				}
				return true;
			}
		}
	}

	public String[] getFragpipeEnzyme() {
		return fragpipeEnzyme;
	}

	public void setFragpipeEnzyme(String[] fragpipeEnzyme) {
		this.fragpipeEnzyme = fragpipeEnzyme;
	}

	public MagDbItem[][] getMagDbItems() {
		return magDbItems;
	}

	public void setMagDbItems(MagDbItem[][] magDbItems) {
		this.magDbItems = magDbItems;
	}

	public MagDbItem getUsedMagDbItem() {
		return usedMagDbItem;
	}

	public void setUsedMagDbItem(MagDbItem usedMagDbItem) {
		this.usedMagDbItem = usedMagDbItem;
	}

	public String getMagDb() {
		return magDb;
	}

	public String getMagDbVersion() {
		return magDbVersion;
	}

	public MagDbItem[] getAvailableMagDbItem() {
		ArrayList<MagDbItem> list = new ArrayList<MagDbItem>();
		if (magDbItems != null && magDbItems.length > 0) {
			for (int i = 0; i < magDbItems.length; i++) {
				for (int j = 0; j < magDbItems[i].length; j++) {
					boolean use = false;
					ArrayList<MagVersion> verlist = magDbItems[i][j].getVersionList();
					for (int k = 0; k < verlist.size(); k++) {
						if (verlist.get(k).isLocalAvailable()) {
							use = true;
							break;
						}
					}

					if (use) {
						list.add(magDbItems[i][j]);
					}
				}
			}
		}
		MagDbItem[] items = list.toArray(new MagDbItem[list.size()]);
		return items;
	}

	public MaxquantModification[] getFullFixMods() {
		return fullFixMods;
	}

	public void setFullFixMods(MaxquantModification[] fullFixMods) {
		this.fullFixMods = fullFixMods;
	}

	public MaxquantModification[] getFullVariMods() {
		return fullVariMods;
	}

	public void setFullVariMods(MaxquantModification[] fullVariMods) {
		this.fullVariMods = fullVariMods;
	}

	// Getters and setters for inherited fields from MetaParameterPFind
	public String[] getFixMods() {
		return fixMods;
	}

	public void setFixMods(String[] fixMods) {
		this.fixMods = fixMods;
	}

	public String[] getVariMods() {
		return variMods;
	}

	public void setVariMods(String[] variMods) {
		this.variMods = variMods;
	}

	public String getEnzyme() {
		return enzyme;
	}

	public void setEnzyme(String enzyme) {
		this.enzyme = enzyme;
	}

	public int getMissCleavages() {
		return missCleavages;
	}

	public void setMissCleavages(int missCleavages) {
		this.missCleavages = missCleavages;
	}

	public int getDigestMode() {
		return digestMode;
	}

	public void setDigestMode(int digestMode) {
		this.digestMode = digestMode;
	}

	public String getQuanMode() {
		return quanMode;
	}

	public void setQuanMode(String quanMode) {
		this.quanMode = quanMode;
	}

	public MaxquantModification[] getIsobaric() {
		return isobaric;
	}

	public void setIsobaric(MaxquantModification[] isobaric) {
		this.isobaric = isobaric;
	}

	public IsobaricTag getIsobaricTag() {
		return isobaricTag;
	}

	public void setIsobaricTag(IsobaricTag isobaricTag) {
		this.isobaricTag = isobaricTag;
	}

	public boolean isCombineLabel() {
		return combineLabel;
	}

	public void setCombineLabel(boolean combineLabel) {
		this.combineLabel = combineLabel;
	}

	public double[][] getIsoCorFactor() {
		return isoCorFactor;
	}

	public void setIsoCorFactor(double[][] isoCorFactor) {
		this.isoCorFactor = isoCorFactor;
	}

	public MetaOpenPar getOsp() {
		return osp;
	}

	public void setOsp(MetaOpenPar osp) {
		this.osp = osp;
	}

	public MetaPep2TaxaPar getPtp() {
		return ptp;
	}

	public void setPtp(MetaPep2TaxaPar ptp) {
		this.ptp = ptp;
	}

	// Methods that delegate to osp (MetaOpenPar)
	public MetaOpenPar getMop() {
		return osp;
	}

	public void setMop(MetaOpenPar mop) {
		this.osp = mop;
	}

	public boolean isOpenSearch() {
		return osp != null ? osp.isOpenSearch() : false;
	}

}
