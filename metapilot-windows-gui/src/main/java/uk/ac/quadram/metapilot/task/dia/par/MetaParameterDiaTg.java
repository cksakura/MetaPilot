package uk.ac.quadram.metapilot.task.dia.par;

import uk.ac.quadram.metapilot.task.par.MetaBaseParameter;
import uk.ac.quadram.metapilot.task.mag.MagDbItem;

public class MetaParameterDiaTg extends MetaParameterDia {

	// Additional fields for DIA Target workflow
	private String[] selectedTaxa;
	private String[] selectedFunctions;
	private String targetLibraryName;
	private String taxonomicLevel;
	private String functionalLevel;
	private String libraryPath;

	public MetaParameterDiaTg() {
		super();
		this.setWorkflowType(uk.ac.quadram.metapilot.MetaPilotApp.WorkflowType.DIA_TARGET.getAnalysisType());
	}
	
	public MetaParameterDiaTg(MetaBaseParameter metaPar, MagDbItem magDbItem, boolean librarySearch,
			boolean selfModelSearch, String[] library, boolean useModel, String[] models) {
		super(metaPar, magDbItem, librarySearch, selfModelSearch, library, useModel, models);
		
		// Initialize additional fields with default values
		this.selectedTaxa = new String[0];
		this.selectedFunctions = new String[0];
		this.targetLibraryName = "target_library";
		this.setWorkflowType(uk.ac.quadram.metapilot.MetaPilotApp.WorkflowType.DIA_TARGET.getAnalysisType());
	}

	public String getLibraryPath() {
		return libraryPath;
	}

	public void setLibraryPath(String libraryPath) {
		this.libraryPath = libraryPath;
	}

	// Getters and setters for additional fields
	public String[] getSelectedTaxa() {
		return selectedTaxa;
	}

	public void setSelectedTaxa(String[] selectedTaxa) {
		this.selectedTaxa = selectedTaxa;
	}

	public String[] getSelectedFunctions() {
		return selectedFunctions;
	}

	public void setSelectedFunctions(String[] selectedFunctions) {
		this.selectedFunctions = selectedFunctions;
	}

	public String getTargetLibraryName() {
		return targetLibraryName;
	}

	public void setTargetLibraryName(String targetLibraryName) {
		this.targetLibraryName = targetLibraryName;
	}

	public String getTaxonomicLevel() {
		return taxonomicLevel;
	}

	public void setTaxonomicLevel(String taxonomicLevel) {
		this.taxonomicLevel = taxonomicLevel;
	}

	public String getFunctionalLevel() {
		return functionalLevel;
	}

	public void setFunctionalLevel(String functionalLevel) {
		this.functionalLevel = functionalLevel;
	}
}
