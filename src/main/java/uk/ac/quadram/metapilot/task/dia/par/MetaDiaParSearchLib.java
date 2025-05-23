package uk.ac.quadram.metapilot.task.dia.par;

import uk.ac.quadram.metapilot.dbSearch.diann.DiannParameter;
import uk.ac.quadram.metapilot.task.dia.DiaLibSearchPar;
import uk.ac.quadram.metapilot.task.mag.par.MetaParameterMag;
import uk.ac.quadram.metapilot.task.par.MetaBaseParameter;

public class MetaDiaParSearchLib extends MetaParameterMag {

	private DiaLibSearchPar diaLibSearchPar;
	private MetaDiaParGenLib diaParGenLib;
	private boolean isTimsTof;

	public MetaDiaParSearchLib(MetaBaseParameter metaPar, String magDb, String magDbVersion) {
		super(metaPar);
		// TODO Auto-generated constructor stub
		this.magDb = magDb;
		this.magDbVersion = magDbVersion;
		this.diaLibSearchPar = new DiaLibSearchPar();
	}

	public MetaDiaParSearchLib(MetaBaseParameter metaPar, String magDb, String magDbVersion,
			DiaLibSearchPar diaLibSearchPar) {
		super(metaPar);
		// TODO Auto-generated constructor stub
		this.magDb = magDb;
		this.magDbVersion = magDbVersion;
		this.diaLibSearchPar = diaLibSearchPar;
	}

	public MetaDiaParSearchLib(MetaBaseParameter metaPar, String magDb, String magDbVersion,
			DiaLibSearchPar diaLibSearchPar, MetaDiaParGenLib diaParGenLib) {
		super(metaPar);
		// TODO Auto-generated constructor stub
		this.magDb = magDb;
		this.magDbVersion = magDbVersion;
		this.diaLibSearchPar = diaLibSearchPar;
		this.diaParGenLib = diaParGenLib;
	}

	public DiaLibSearchPar getDiaLibSearchPar() {
		return diaLibSearchPar;
	}

	public void setDiaLibSearchPar(DiaLibSearchPar diaLibSearchPar) {
		this.diaLibSearchPar = diaLibSearchPar;
	}

	public DiannParameter getDiannPar() {
		return diaParGenLib.getDiannPar();
	}

	public MetaDiaParGenLib getDiaParGenLib() {
		return diaParGenLib;
	}

	public void setDiaParGenLib(MetaDiaParGenLib diaParGenLib) {
		this.diaParGenLib = diaParGenLib;
	}

	public boolean isTimsTof() {
		return isTimsTof;
	}

	public void setTimsTof(boolean isTimsTof) {
		this.isTimsTof = isTimsTof;
	}

}
