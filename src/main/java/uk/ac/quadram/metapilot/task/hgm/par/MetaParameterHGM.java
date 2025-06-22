package uk.ac.quadram.metapilot.task.hgm.par;

import uk.ac.quadram.metapilot.core.mod.IsobaricTag;
import uk.ac.quadram.metapilot.dbSearch.MaxQuant.MaxquantModification;
import uk.ac.quadram.metapilot.task.par.MetaBaseParameter;
import uk.ac.quadram.metapilot.task.par.MetaPep2TaxaPar;
import uk.ac.quadram.metapilot.task.pfind.par.MetaParameterPFind;
import uk.ac.quadram.metapilot.task.v2.par.MetaOpenPar;

public class MetaParameterHGM extends MetaParameterPFind {

	public MetaParameterHGM(MetaBaseParameter metaPar) {
		super(metaPar);
	}

	public MetaParameterHGM(MetaBaseParameter metaPar, String[] fixMods, String[] variMods, String enzyme,
			int missCleavages, int digestMode, String quanMode, MaxquantModification[] isobaric,
			IsobaricTag isobaricTag, boolean combineLabel, double[][] isoCorFactor) {
		super(metaPar, fixMods, variMods, enzyme, missCleavages, digestMode, quanMode, isobaric, isobaricTag,
				combineLabel, isoCorFactor);
		this.osp = MetaOpenPar.getDefault();
		this.ptp = MetaPep2TaxaPar.getDefault();
	}

}