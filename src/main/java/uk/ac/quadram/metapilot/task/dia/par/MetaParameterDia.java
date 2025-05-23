package uk.ac.quadram.metapilot.task.dia.par;

import uk.ac.quadram.metapilot.task.mag.par.MetaParameterMag;
import uk.ac.quadram.metapilot.task.par.MetaBaseParameter;

public class MetaParameterDia extends MetaParameterMag {

	private boolean librarySearch;
	private boolean selfModelSearch;
	private boolean useModel;
	private String[] library;
	private String[] models;

	public MetaParameterDia(MetaBaseParameter metaPar, String magDb, String magDbVersion, boolean librarySearch,
			boolean selfModelSearch, String[] library, boolean useModel, String[] models) {
		super(metaPar);
		// TODO Auto-generated constructor stub
		this.magDb = magDb;
		this.magDbVersion = magDbVersion;
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

}
