package uk.ac.quadram.metapilot.core.runner;

import uk.ac.quadram.metapilot.core.function.v2.EnzymeCommission;
import uk.ac.quadram.metapilot.core.function.v2.GoObo;
import uk.ac.quadram.metapilot.task.mag.MetaProteinAnnoMag;

import java.util.HashMap;

/**
 * In-memory MAG functional annotation result.
 */
public final class FunctionalAnnotationResult {
	public MetaProteinAnnoMag[] annotatedProteins = new MetaProteinAnnoMag[0];
	public HashMap<String, String[]> usedCogMap = new HashMap<>();
	public HashMap<String, String[]> usedNogMap = new HashMap<>();
	public HashMap<String, String> usedKeggMap = new HashMap<>();
	public HashMap<String, GoObo> usedGoMap = new HashMap<>();
	public HashMap<String, EnzymeCommission> usedEcMap = new HashMap<>();
}
