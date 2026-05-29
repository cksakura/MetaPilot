package uk.ac.quadram.metapilot.dbSearch.sage;

import java.io.BufferedInputStream;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Date;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;

import javax.swing.JProgressBar;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.json.JSONArray;
import org.json.JSONObject;

import uk.ac.quadram.metapilot.task.dda.igc.par.DDAIgcSageParameters;
import uk.ac.quadram.metapilot.task.dda.mag.par.DDAMagSageParameters;

public class SageTask {

	private static String taskName = "Sage task";
	private static final Logger LOGGER = LogManager.getLogger(SageTask.class);
	private SimpleDateFormat format = new SimpleDateFormat("MM_dd_yyyy_HH_mm_ss");

	private File sageFile;
	private ArrayList<String> parList;
	private StringBuilder headsb;
	private StringBuilder tailsb;

	public SageTask(String sage) {
		this(new File(sage));
	}

	public SageTask(File sage) {
		this.sageFile = sage;
		this.parList = new ArrayList<String>();
		this.headsb = new StringBuilder();
		this.headsb.append("@echo off\n");
		this.headsb.append("cd /d \"" + sageFile.getParent() + "\"\n");

		this.tailsb = new StringBuilder();
		this.tailsb.append("if errorlevel 1 (\n");
		this.tailsb.append("\techo Error: Failed to start the application.\n");
		this.tailsb.append(") else (\n");
		this.tailsb.append("\techo Application started successfully.\n");
		this.tailsb.append(")\n");
		this.tailsb.append("exit");
	}

	public void addTask(File configFile) {
		File batFile = new File(configFile.getParentFile(), "sage.bat");
		try (PrintWriter writer = new PrintWriter(batFile)) {
			writer.println(headsb);
			writer.println("\"" + sageFile.getAbsolutePath() + "\" \"" + configFile.getAbsolutePath() + "\"");
			writer.println();
			writer.print("exit");
			writer.close();
		} catch (IOException e) {
			LOGGER.error(taskName + ": error in writing the .bat file to " + batFile, e);
			System.out.println(
					format.format(new Date()) + "\t" + taskName + ": error in writing the .bat file to " + batFile);
		}
		this.parList.add(batFile.getAbsolutePath());
	}

	/**
	 * Search mzML files against a FASTA database using Sage
	 * 
	 * @param mzMLFiles       Array of mzML file paths
	 * @param result          Output directory for results
	 * @param db              FASTA database file path
	 * @param precursorTol    Precursor tolerance in ppm
	 * @param fragmentTol     Fragment tolerance in ppm
	 * @param enzyme          Enzyme name
	 * @param missedCleavages Number of missed cleavages
	 * @param minPepLen       Minimum peptide length
	 * @param maxPepLen       Maximum peptide length
	 * @param fdr             FDR threshold
	 * @param fixedMods       Array of fixed modifications
	 * @param variableMods    Array of variable modifications
	 */
	public void searchMzML(String[] mzMLFiles, String result, String db, double precursorTol, double fragmentTol,
			String enzyme, int missedCleavages, int minPepLen, int maxPepLen, double fdr, String[] fixedMods,
			String[] variableMods) {
		
		File resultDir = new File(result);
		if (!resultDir.exists()) {
			resultDir.mkdirs();
		}

		File configFile = new File(resultDir, "sage_config.json");
		
		try (PrintWriter writer = new PrintWriter(configFile)) {
			// Create JSON structure using JSONObject
			JSONObject config = new JSONObject();
			
			// Database section
			JSONObject database = new JSONObject();
			database.put("bucket_size", 8192);
			
			// Enzyme section
			JSONObject enzymeObj = new JSONObject();
			enzymeObj.put("missed_cleavages", missedCleavages);
			enzymeObj.put("min_len", minPepLen);
			enzymeObj.put("max_len", maxPepLen);
			
			// Only add cleave_at if it's not empty
			String cleaveAt = getEnzymeCleavageSites(enzyme);
			if (!cleaveAt.isEmpty()) {
				enzymeObj.put("cleave_at", cleaveAt);
			}
			
			// Only add restrict if it's not empty
			String restrict = getEnzymeRestriction(enzyme);
			if (!restrict.isEmpty()) {
				enzymeObj.put("restrict", restrict);
			}
			
			database.put("enzyme", enzymeObj);
			database.put("peptide_min_mass", 500.0);
			database.put("peptide_max_mass", 5000.0);
			database.put("ion_kinds", new JSONArray().put("b").put("y"));
			database.put("min_ion_index", 2);
			database.put("max_variable_mods", 3);
			
			// Fixed modifications
			JSONObject staticMods = new JSONObject();
			if (fixedMods != null && fixedMods.length > 0) {
				for (String mod : fixedMods) {
					if (mod != null && !mod.trim().isEmpty()) {
						String[] parts = mod.split("\\[|\\]");
						if (parts.length >= 2) {
							String aminoAcid = parts[1]; // Extract amino acid from [C]
							String massShift = getModificationMass(parts[0]); // Get mass for modification name
							staticMods.put(aminoAcid, Double.parseDouble(massShift));
						}
					}
				}
			}
			database.put("static_mods", staticMods);
			
			// Variable modifications
			JSONObject variableModsObj = new JSONObject();
			if (variableMods != null && variableMods.length > 0) {
				for (String mod : variableMods) {
					if (mod != null && !mod.trim().isEmpty()) {
						String[] parts = mod.split("\\[|\\]");
						if (parts.length >= 2) {
							String aminoAcid = parts[1]; // Extract amino acid from [M]
							String massShift = getModificationMass(parts[0]); // Get mass for modification name
							variableModsObj.put(aminoAcid, Double.parseDouble(massShift));
						}
					}
				}
			}
			database.put("variable_mods", variableModsObj);
			database.put("decoy_tag", "REV_");
			database.put("generate_decoys", true);
			database.put("fasta", db);
			
			config.put("database", database);
			
			// Precursor tolerance (open search format - use da)
			JSONObject precursorTolObj = new JSONObject();
			precursorTolObj.put("da", new JSONArray().put(-precursorTol).put(precursorTol));
			config.put("precursor_tol", precursorTolObj);
			
			// Fragment tolerance
			JSONObject fragmentTolObj = new JSONObject();
			fragmentTolObj.put("ppm", new JSONArray().put(-fragmentTol).put(fragmentTol));
			config.put("fragment_tol", fragmentTolObj);
			
			// Open search specific parameters
			config.put("deisotope", true);
			config.put("chimera", true);
			config.put("min_peaks", 15);
			config.put("max_peaks", 250);
			config.put("min_matched_peaks", 4);
			config.put("report_psms", 2);
			config.put("predict_rt", true);
			
			// Output directory
			config.put("output_directory", result);
			
			// mzML paths
			JSONArray mzmlPaths = new JSONArray();
			for (String mzmlFile : mzMLFiles) {
				mzmlPaths.put(mzmlFile);
			}
			config.put("mzml_paths", mzmlPaths);
			
			// Write JSON to file
			writer.println(config.toString(2)); // Pretty print with 2-space indentation

		} catch (IOException e) {
			LOGGER.error(taskName + ": error in writing Sage config file to " + configFile, e);
			System.out.println(format.format(new Date()) + "\t" + taskName + ": error in writing Sage config file to "
					+ configFile);
			return;
		}

		addTask(configFile);
	}
	
	/**
	 * Search mzML files against a FASTA database using Sage with enhanced
	 * parameters from DDAMagSageParameters object
	 * 
	 * @param mzMLFiles       Array of mzML file paths
	 * @param result          Output directory for results
	 * @param db              FASTA database file path
	 * @param parameter       DDAMagSageParameters object containing all parameters
	 * @param isQuant         Whether to include quantification settings
	 */
	public void searchMzML(String[] mzMLFiles, String result, String db, DDAIgcSageParameters parameter,
			boolean isQuant) {

		// Extract all values from the parameter object
		int missedCleavages = parameter.getMissCleavages();
		String enzyme = parameter.getEnzyme();
		String[] fixedMods = parameter.getFixMods();
		String[] variableMods = parameter.getVariMods();
		double precursorTolMin = parameter.getPrecursorToleranceMin();
		double precursorTolMax = parameter.getPrecursorToleranceMax();
		double fragmentTol = parameter.getFragmentTolerance();
		String precursorUnit = parameter.getPrecursorToleranceUnit();
		String fragmentUnit = parameter.getFragmentToleranceUnit();
		int minPepLen = parameter.getMinPeptideLength();
		int maxPepLen = parameter.getMaxPeptideLength();

		File resultDir = new File(result);
		if (!resultDir.exists()) {
			resultDir.mkdirs();
		}

		File configFile = new File(resultDir, "sage_config.json");

		try (PrintWriter writer = new PrintWriter(configFile)) {
			// Create JSON structure using JSONObject
			JSONObject config = new JSONObject();
			
			// Database section
			JSONObject database = new JSONObject();
			database.put("bucket_size", 8192);
			
			// Enzyme section
			JSONObject enzymeObj = new JSONObject();
			enzymeObj.put("missed_cleavages", missedCleavages);
			enzymeObj.put("min_len", minPepLen);
			enzymeObj.put("max_len", maxPepLen);
			
			// Only add cleave_at if it's not empty
			String cleaveAt = getEnzymeCleavageSites(enzyme);
			if (!cleaveAt.isEmpty()) {
				enzymeObj.put("cleave_at", cleaveAt);
			}
			
			// Only add restrict if it's not empty
			String restrict = getEnzymeRestriction(enzyme);
			if (!restrict.isEmpty()) {
				enzymeObj.put("restrict", restrict);
			}
			
			database.put("enzyme", enzymeObj);
			database.put("peptide_min_mass", 500.0);
			database.put("peptide_max_mass", 5000.0);
			database.put("ion_kinds", new JSONArray().put("b").put("y"));
			database.put("min_ion_index", 2);
			database.put("max_variable_mods", 3);
			
			// Fixed modifications - Format: "site:mass" (e.g., "C:57.02146", "peptide n-term:42.0106")
			JSONObject staticMods = new JSONObject();
			if (fixedMods != null && fixedMods.length > 0) {
				for (String mod : fixedMods) {
					if (mod != null && !mod.trim().isEmpty()) {
						String[] parts = mod.split(":");
						if (parts.length >= 2) {
							String site = parts[0].trim(); // Extract site (e.g., "C", "peptide n-term")
							String massShift = parts[1].trim(); // Extract mass value
							
							// Convert site to Sage notation (handles terminal modifications)
							String sageKey = convertSiteToSageNotation(site);
							
							staticMods.put(sageKey, Double.parseDouble(massShift));
						}
					}
				}
			}
			database.put("static_mods", staticMods);
			
			// Variable modifications - Format: "modName:sites:mass" (e.g., "Oxidation:M:15.9949", "Acetylation:protein n-term:42.0106")
			JSONObject variableModsObj = new JSONObject();
			if (variableMods != null && variableMods.length > 0) {
				for (String mod : variableMods) {
					if (mod != null && !mod.trim().isEmpty()) {
						String[] parts = mod.split(":");
						if (parts.length >= 3) {
							// parts[0] is modName (not used in Sage config)
							String sites = parts[1].trim(); // Sites (e.g., "M", "STY", "protein n-term")
							String massShift = parts[2].trim(); // Extract mass value
							
							// Handle multiple sites (e.g., "STY" -> "S", "T", "Y")
							if (sites.length() > 1 && !sites.contains(" ") && !sites.contains("-")) {
								// Multiple amino acids like "STY"
								for (char aminoAcid : sites.toCharArray()) {
									String sageKey = convertSiteToSageNotation(String.valueOf(aminoAcid));
									// Sage variable mods use array format
									if (variableModsObj.has(sageKey)) {
										// Append to existing array
										JSONArray existingArray = variableModsObj.getJSONArray(sageKey);
										existingArray.put(Double.parseDouble(massShift));
									} else {
										// Create new array
										variableModsObj.put(sageKey, new JSONArray().put(Double.parseDouble(massShift)));
									}
								}
							} else {
								// Single site or terminal modification
								String sageKey = convertSiteToSageNotation(sites);
								
								// Sage variable mods use array format
								if (variableModsObj.has(sageKey)) {
									// Append to existing array
									JSONArray existingArray = variableModsObj.getJSONArray(sageKey);
									existingArray.put(Double.parseDouble(massShift));
								} else {
									// Create new array
									variableModsObj.put(sageKey, new JSONArray().put(Double.parseDouble(massShift)));
								}
							}
						}
					}
				}
			}
			database.put("variable_mods", variableModsObj);
			database.put("decoy_tag", "REV_");
			database.put("generate_decoys", true);
			database.put("fasta", db);
			
			config.put("database", database);
			
			if (isQuant) {
				// Quant format - add quant section first
				JSONObject quant = new JSONObject();
				quant.put("lfq", true);
				
				JSONObject lfqSettings = new JSONObject();
				lfqSettings.put("peak_scoring", parameter.getPeakScoring());
				lfqSettings.put("integration", parameter.getIntegration());
				lfqSettings.put("spectral_angle", parameter.getSpectralAngle());
				lfqSettings.put("ppm_tolerance", parameter.getLfqPpmTolerance());
				quant.put("lfq_settings", lfqSettings);
				
				config.put("quant", quant);
				
				// Quant format - use unit from parameter (usually "ppm" for both)
				JSONObject precursorTolObj = new JSONObject();
				precursorTolObj.put(precursorUnit.toLowerCase(), new JSONArray().put(precursorTolMin).put(precursorTolMax));
				config.put("precursor_tol", precursorTolObj);
				
				JSONObject fragmentTolObj = new JSONObject();
				fragmentTolObj.put(fragmentUnit.toLowerCase(), new JSONArray().put(-fragmentTol).put(fragmentTol));
				config.put("fragment_tol", fragmentTolObj);
				
				config.put("isotope_errors", new JSONArray().put(0).put(2));
				config.put("deisotope", true);
				config.put("min_peaks", 15);
				config.put("max_peaks", 150);
				config.put("max_fragment_charge", 1);
				config.put("min_matched_peaks", 4);
				config.put("predict_rt", true);
			} else {
				// Open search format - use unit from parameter (usually "da" for precursor, "ppm" for fragment)
				JSONObject precursorTolObj = new JSONObject();
				precursorTolObj.put(precursorUnit.toLowerCase(), new JSONArray().put(precursorTolMin).put(precursorTolMax));
				config.put("precursor_tol", precursorTolObj);
				
				JSONObject fragmentTolObj = new JSONObject();
				fragmentTolObj.put(fragmentUnit.toLowerCase(), new JSONArray().put(-fragmentTol).put(fragmentTol));
				config.put("fragment_tol", fragmentTolObj);
				
				config.put("deisotope", true);
				config.put("chimera", true);
				config.put("min_peaks", 15);
				config.put("max_peaks", 250);
				config.put("min_matched_peaks", 4);
				config.put("report_psms", 2);
				config.put("predict_rt", true);
			}
			
			// Output directory (same for both formats)
			config.put("output_directory", result);
			
			// mzML paths (same for both formats)
			JSONArray mzmlPaths = new JSONArray();
			for (String mzmlFile : mzMLFiles) {
				mzmlPaths.put(mzmlFile);
			}
			config.put("mzml_paths", mzmlPaths);
			
			// Write JSON to file
			writer.println(config.toString(2)); // Pretty print with 2-space indentation
			
		} catch (IOException e) {
			LOGGER.error(taskName + ": error in writing Sage config file to " + configFile, e);
			System.out.println(format.format(new Date()) + "\t" + taskName + ": error in writing Sage config file to "
					+ configFile);
			return;
		}
		
		addTask(configFile);
	}

	/**
	 * Search mzML files against a FASTA database using Sage with enhanced
	 * parameters from DDAMagSageParameters object
	 * 
	 * @param mzMLFiles       Array of mzML file paths
	 * @param result          Output directory for results
	 * @param db              FASTA database file path
	 * @param parameter       DDAMagSageParameters object containing all parameters
	 * @param isQuant         Whether to include quantification settings
	 */
	public void searchMzML(String[] mzMLFiles, String result, String db, DDAMagSageParameters parameter,
			boolean isQuant) {

		// Extract all values from the parameter object
		int missedCleavages = parameter.getMissCleavages();
		String enzyme = parameter.getEnzyme();
		String[] fixedMods = parameter.getFixMods();
		String[] variableMods = parameter.getVariMods();
		double precursorTolMin = parameter.getPrecursorToleranceMin();
		double precursorTolMax = parameter.getPrecursorToleranceMax();
		double fragmentTol = parameter.getFragmentTolerance();
		String precursorUnit = parameter.getPrecursorToleranceUnit();
		String fragmentUnit = parameter.getFragmentToleranceUnit();
		int minPepLen = parameter.getMinPeptideLength();
		int maxPepLen = parameter.getMaxPeptideLength();

		File resultDir = new File(result);
		if (!resultDir.exists()) {
			resultDir.mkdirs();
		}

		File configFile = new File(resultDir, "sage_config.json");

		try (PrintWriter writer = new PrintWriter(configFile)) {
			// Create JSON structure using JSONObject
			JSONObject config = new JSONObject();
			
			// Database section
			JSONObject database = new JSONObject();
			database.put("bucket_size", 8192);
			
			// Enzyme section
			JSONObject enzymeObj = new JSONObject();
			enzymeObj.put("missed_cleavages", missedCleavages);
			enzymeObj.put("min_len", minPepLen);
			enzymeObj.put("max_len", maxPepLen);
			
			// Only add cleave_at if it's not empty
			String cleaveAt = getEnzymeCleavageSites(enzyme);
			if (!cleaveAt.isEmpty()) {
				enzymeObj.put("cleave_at", cleaveAt);
			}
			
			// Only add restrict if it's not empty
			String restrict = getEnzymeRestriction(enzyme);
			if (!restrict.isEmpty()) {
				enzymeObj.put("restrict", restrict);
			}
			
			database.put("enzyme", enzymeObj);
			database.put("peptide_min_mass", 500.0);
			database.put("peptide_max_mass", 5000.0);
			database.put("ion_kinds", new JSONArray().put("b").put("y"));
			database.put("min_ion_index", 2);
			database.put("max_variable_mods", 3);
			
			// Fixed modifications - Format: "site:mass" (e.g., "C:57.02146", "peptide n-term:42.0106")
			JSONObject staticMods = new JSONObject();
			if (fixedMods != null && fixedMods.length > 0) {
				for (String mod : fixedMods) {
					if (mod != null && !mod.trim().isEmpty()) {
						String[] parts = mod.split(":");
						if (parts.length >= 2) {
							String site = parts[0].trim(); // Extract site (e.g., "C", "peptide n-term")
							String massShift = parts[1].trim(); // Extract mass value
							
							// Convert site to Sage notation (handles terminal modifications)
							String sageKey = convertSiteToSageNotation(site);
							
							staticMods.put(sageKey, Double.parseDouble(massShift));
						}
					}
				}
			}
			database.put("static_mods", staticMods);
			
			// Variable modifications - Format: "modName:sites:mass" (e.g., "Oxidation:M:15.9949", "Acetylation:protein n-term:42.0106")
			JSONObject variableModsObj = new JSONObject();
			if (variableMods != null && variableMods.length > 0) {
				for (String mod : variableMods) {
					if (mod != null && !mod.trim().isEmpty()) {
						String[] parts = mod.split(":");
						if (parts.length >= 3) {
							// parts[0] is modName (not used in Sage config)
							String sites = parts[1].trim(); // Sites (e.g., "M", "STY", "protein n-term")
							String massShift = parts[2].trim(); // Extract mass value
							
							// Handle multiple sites (e.g., "STY" -> "S", "T", "Y")
							if (sites.length() > 1 && !sites.contains(" ") && !sites.contains("-")) {
								// Multiple amino acids like "STY"
								for (char aminoAcid : sites.toCharArray()) {
									String sageKey = convertSiteToSageNotation(String.valueOf(aminoAcid));
									// Sage variable mods use array format
									if (variableModsObj.has(sageKey)) {
										// Append to existing array
										JSONArray existingArray = variableModsObj.getJSONArray(sageKey);
										existingArray.put(Double.parseDouble(massShift));
									} else {
										// Create new array
										variableModsObj.put(sageKey, new JSONArray().put(Double.parseDouble(massShift)));
									}
								}
							} else {
								// Single site or terminal modification
								String sageKey = convertSiteToSageNotation(sites);
								
								// Sage variable mods use array format
								if (variableModsObj.has(sageKey)) {
									// Append to existing array
									JSONArray existingArray = variableModsObj.getJSONArray(sageKey);
									existingArray.put(Double.parseDouble(massShift));
								} else {
									// Create new array
									variableModsObj.put(sageKey, new JSONArray().put(Double.parseDouble(massShift)));
								}
							}
						}
					}
				}
			}
			database.put("variable_mods", variableModsObj);
			database.put("decoy_tag", "REV_");
			database.put("generate_decoys", true);
			database.put("fasta", db);
			
			config.put("database", database);
			
			if (isQuant) {
				// Quant format - add quant section first
				JSONObject quant = new JSONObject();
				quant.put("lfq", true);
				
				JSONObject lfqSettings = new JSONObject();
				lfqSettings.put("peak_scoring", parameter.getPeakScoring());
				lfqSettings.put("integration", parameter.getIntegration());
				lfqSettings.put("spectral_angle", parameter.getSpectralAngle());
				lfqSettings.put("ppm_tolerance", parameter.getLfqPpmTolerance());
				quant.put("lfq_settings", lfqSettings);
				
				config.put("quant", quant);
				
				// Quant format - use unit from parameter (usually "ppm" for both)
				JSONObject precursorTolObj = new JSONObject();
				precursorTolObj.put(precursorUnit.toLowerCase(), new JSONArray().put(precursorTolMin).put(precursorTolMax));
				config.put("precursor_tol", precursorTolObj);
				
				JSONObject fragmentTolObj = new JSONObject();
				fragmentTolObj.put(fragmentUnit.toLowerCase(), new JSONArray().put(-fragmentTol).put(fragmentTol));
				config.put("fragment_tol", fragmentTolObj);
				
				config.put("isotope_errors", new JSONArray().put(0).put(2));
				config.put("deisotope", true);
				config.put("min_peaks", 15);
				config.put("max_peaks", 150);
				config.put("max_fragment_charge", 1);
				config.put("min_matched_peaks", 4);
				config.put("predict_rt", true);
			} else {
				// Open search format - use unit from parameter (usually "da" for precursor, "ppm" for fragment)
				JSONObject precursorTolObj = new JSONObject();
				precursorTolObj.put(precursorUnit.toLowerCase(), new JSONArray().put(precursorTolMin).put(precursorTolMax));
				config.put("precursor_tol", precursorTolObj);
				
				JSONObject fragmentTolObj = new JSONObject();
				fragmentTolObj.put(fragmentUnit.toLowerCase(), new JSONArray().put(-fragmentTol).put(fragmentTol));
				config.put("fragment_tol", fragmentTolObj);
				
				config.put("deisotope", true);
				config.put("chimera", true);
				config.put("min_peaks", 15);
				config.put("max_peaks", 250);
				config.put("min_matched_peaks", 4);
				config.put("report_psms", 2);
				config.put("predict_rt", true);
			}
			
			// Output directory (same for both formats)
			config.put("output_directory", result);
			
			// mzML paths (same for both formats)
			JSONArray mzmlPaths = new JSONArray();
			for (String mzmlFile : mzMLFiles) {
				mzmlPaths.put(mzmlFile);
			}
			config.put("mzml_paths", mzmlPaths);
			
			// Write JSON to file
			writer.println(config.toString(2)); // Pretty print with 2-space indentation

		} catch (IOException e) {
			LOGGER.error(taskName + ": error in writing Sage config file to " + configFile, e);
			System.out.println(format.format(new Date()) + "\t" + taskName + ": error in writing Sage config file to "
					+ configFile);
			return;
		}

		addTask(configFile);
	}
	
	/**
	 * Search mzML files against a FASTA database using Sage with enhanced
	 * parameters
	 * 
	 * @param mzMLFiles       Array of mzML file paths
	 * @param result          Output directory for results
	 * @param db              FASTA database file path
	 * @param precursorTolMin Minimum precursor tolerance in ppm
	 * @param precursorTolMax Maximum precursor tolerance in ppm
	 * @param fragmentTol     Fragment tolerance in ppm
	 * @param enzyme          Enzyme name
	 * @param missedCleavages Number of missed cleavages
	 * @param fixedMods       Array of fixed modifications
	 * @param variableMods    Array of variable modifications
	 * @param isQuant         Whether to include quantification settings
	 */
	public void searchMzML(String[] mzMLFiles, String result, String db, double precursorTolMin, double precursorTolMax,
			double fragmentTol, String enzyme, int missedCleavages, String[] fixedMods, String[] variableMods,
			boolean isQuant) {
		
		File resultDir = new File(result);
		if (!resultDir.exists()) {
			resultDir.mkdirs();
		}

		File configFile = new File(resultDir, "sage_config.json");
		
		try (PrintWriter writer = new PrintWriter(configFile)) {
			// Create JSON structure using JSONObject
			JSONObject config = new JSONObject();
			
			// Database section
			JSONObject database = new JSONObject();
			database.put("bucket_size", 8192);
			
			// Enzyme section
			JSONObject enzymeObj = new JSONObject();
			enzymeObj.put("missed_cleavages", missedCleavages);
			enzymeObj.put("min_len", 7);
			enzymeObj.put("max_len", 50);
			
			// Only add cleave_at if it's not empty
			String cleaveAt = getEnzymeCleavageSites(enzyme);
			if (!cleaveAt.isEmpty()) {
				enzymeObj.put("cleave_at", cleaveAt);
			}
			
			// Only add restrict if it's not empty
			String restrict = getEnzymeRestriction(enzyme);
			if (!restrict.isEmpty()) {
				enzymeObj.put("restrict", restrict);
			}
			
			database.put("enzyme", enzymeObj);
			database.put("peptide_min_mass", 500.0);
			database.put("peptide_max_mass", 5000.0);
			database.put("ion_kinds", new JSONArray().put("b").put("y"));
			database.put("min_ion_index", 2);
			database.put("max_variable_mods", 3);
			
			// Fixed modifications
			JSONObject staticMods = new JSONObject();
			if (fixedMods != null && fixedMods.length > 0) {
				for (String mod : fixedMods) {
					if (mod != null && !mod.trim().isEmpty()) {
						String[] parts = mod.split(":");
						if (parts.length >= 2) {
							String aminoAcid = parts[0]; // Extract amino acid
							String massShift = parts[1]; // Extract mass value
							staticMods.put(aminoAcid, Double.parseDouble(massShift));
						}
					}
				}
			}
			database.put("static_mods", staticMods);
			
			// Variable modifications
			JSONObject variableModsObj = new JSONObject();
			if (variableMods != null && variableMods.length > 0) {
				for (String mod : variableMods) {
					if (mod != null && !mod.trim().isEmpty()) {
						String[] parts = mod.split(":");
						if (parts.length >= 3) {
							String aminoAcid = parts[1]; // Extract amino acid from middle part
							String massShift = parts[2]; // Extract mass value from last part
							variableModsObj.put(aminoAcid, Double.parseDouble(massShift));
						}
					}
				}
			}
			database.put("variable_mods", variableModsObj);
			database.put("decoy_tag", "REV_");
			database.put("generate_decoys", true);
			database.put("fasta", db);
			
			config.put("database", database);
			
			if (isQuant) {
				// Quant format - add quant section first
				JSONObject quant = new JSONObject();
				quant.put("lfq", true);
				
				JSONObject lfqSettings = new JSONObject();
				lfqSettings.put("peak_scoring", "Hybrid");
				lfqSettings.put("integration", "Sum");
				lfqSettings.put("spectral_angle", 0.6);
				lfqSettings.put("ppm_tolerance", 5.0);
				quant.put("lfq_settings", lfqSettings);
				
				config.put("quant", quant);
				
				// Quant format - use ppm for both tolerances
				JSONObject precursorTolObj = new JSONObject();
				precursorTolObj.put("ppm", new JSONArray().put(precursorTolMin).put(precursorTolMax));
				config.put("precursor_tol", precursorTolObj);
				
				JSONObject fragmentTolObj = new JSONObject();
				fragmentTolObj.put("ppm", new JSONArray().put(-fragmentTol).put(fragmentTol));
				config.put("fragment_tol", fragmentTolObj);
				
				config.put("isotope_errors", new JSONArray().put(0).put(2));
				config.put("deisotope", true);
				config.put("min_peaks", 15);
				config.put("max_peaks", 150);
				config.put("max_fragment_charge", 1);
				config.put("min_matched_peaks", 4);
				config.put("predict_rt", true);
			} else {
				// Open search format - use da for precursor, ppm for fragment
				JSONObject precursorTolObj = new JSONObject();
				precursorTolObj.put("da", new JSONArray().put(precursorTolMin).put(precursorTolMax));
				config.put("precursor_tol", precursorTolObj);
				
				JSONObject fragmentTolObj = new JSONObject();
				fragmentTolObj.put("ppm", new JSONArray().put(-fragmentTol).put(fragmentTol));
				config.put("fragment_tol", fragmentTolObj);
				
				config.put("deisotope", true);
				config.put("chimera", true);
				config.put("min_peaks", 15);
				config.put("max_peaks", 250);
				config.put("min_matched_peaks", 4);
				config.put("report_psms", 2);
				config.put("predict_rt", true);
			}
			
			// Output directory (same for both formats)
			config.put("output_directory", result);
			
			// mzML paths (same for both formats)
			JSONArray mzmlPaths = new JSONArray();
			for (String mzmlFile : mzMLFiles) {
				mzmlPaths.put(mzmlFile);
			}
			config.put("mzml_paths", mzmlPaths);
			
			// Write JSON to file
			writer.println(config.toString(2)); // Pretty print with 2-space indentation
			
		} catch (IOException e) {
			LOGGER.error(taskName + ": error in writing Sage config file to " + configFile, e);
			System.out.println(format.format(new Date()) + "\t" + taskName + ": error in writing Sage config file to "
					+ configFile);
			return;
		}
		
		addTask(configFile);
	}

	/**
	 * Convert UI modification site to Sage terminal notation
	 * 
	 * @param site The site from UI (e.g., "C", "peptide n-term", "protein n-term")
	 * @return Sage notation (e.g., "C", "^", "[")
	 */
	private String convertSiteToSageNotation(String site) {
		if (site == null || site.trim().isEmpty()) {
			return site;
		}
		
		String siteLower = site.toLowerCase().trim();
		
		// Terminal modifications
		if (siteLower.equals("peptide n-term")) {
			return "^"; // Peptide N-terminus (any amino acid)
		} else if (siteLower.equals("peptide c-term")) {
			return "$"; // Peptide C-terminus
		} else if (siteLower.equals("protein n-term")) {
			return "["; // Protein N-terminus
		} else if (siteLower.equals("protein c-term")) {
			return "]"; // Protein C-terminus
		}
		
		// Check for site-specific peptide N-terminal modifications (e.g., "peptide n-term Q")
		// Note: Current UI doesn't support this, but adding for future compatibility
		if (siteLower.startsWith("peptide n-term ")) {
			String aminoAcid = site.substring("peptide n-term ".length()).trim();
			if (aminoAcid.length() == 1) {
				return "^" + aminoAcid.toUpperCase(); // e.g., "^Q" for peptide N-term Q
			}
		}
		
		// Regular amino acid - return as is
		return site;
	}

	/**
	 * Get the mass shift for a modification name
	 */
	private String getModificationMass(String modificationName) {
		switch (modificationName.trim()) {
			case "Carbamidomethyl":
				return "57.0214";
			case "Oxidation":
				return "15.9949";
			case "Acetyl":
				return "42.0106";
			case "Phospho":
				return "79.9663";
			case "Methyl":
				return "14.0157";
			case "Dimethyl":
				return "28.0314";
			case "Trimethyl":
				return "42.0471";
			case "Deamidated":
				return "0.9840";
			case "Pyro-glu":
				return "-17.0265";
			case "Pyro_glu":
				return "-17.0265";
			default:
				LOGGER.warn("Unknown modification: " + modificationName + ", using default mass 0");
				return "0";
		}
	}

	/**
	 * Get enzyme cleavage sites based on enzyme name
	 */
	public String getEnzymeCleavageSites(String enzyme) {
		if (enzyme == null)
			return "KR";

		switch (enzyme.toLowerCase().trim()) {
		case "trypsin":
		case "trypsin/p":
			return "KR";
		case "lys-c":
		case "lysyl endopeptidase":
			return "K";
		case "arg-c":
		case "arginyl endopeptidase":
			return "R";
		case "glu-c":
		case "glutamyl endopeptidase":
			return "E";
		case "asp-n":
		case "aspartyl endopeptidase":
			return "D";
		case "chymotrypsin":
			return "FWY";
		case "elastase":
			return "AGILV";
		case "thermolysin":
			return "AILMFV";
		case "proteinase k":
			return "AEFILMVWY";
		case "pepsin":
			return "FL";
		case "clostripain":
			return "R";
		case "staph protease":
			return "E";
		case "no enzyme":
		case "unspecific":
			return ""; // Return empty string - will be omitted from JSON
		default:
			LOGGER.warn("Unknown enzyme: " + enzyme + ", using default trypsin cleavage sites (KR)");
			return "KR";
		}
	}

	/**
	 * Get enzyme restriction sites based on enzyme name
	 */
	public String getEnzymeRestriction(String enzyme) {
		if (enzyme == null)
			return "P";

		switch (enzyme.toLowerCase().trim()) {
		case "trypsin/p":
			return "P";
		case "trypsin":
		case "lys-c":
		case "lysyl endopeptidase":
		case "arg-c":
		case "arginyl endopeptidase":
		case "glu-c":
		case "glutamyl endopeptidase":
		case "asp-n":
		case "aspartyl endopeptidase":
		case "chymotrypsin":
		case "elastase":
		case "thermolysin":
		case "proteinase k":
		case "pepsin":
		case "clostripain":
		case "staph protease":
			return ""; // Return empty string - will be omitted from JSON
		case "no enzyme":
		case "unspecific":
			return ""; // Return empty string - will be omitted from JSON
		default:
			LOGGER.warn("Unknown enzyme: " + enzyme + ", using default restriction (P)");
			return "P";
		}
	}

	/**
	 * Set total thread count for parallel execution
	 */
	public void setTotalThread(int threadCount) {
		// Sage handles threading internally, but we can use this for task scheduling
		LOGGER.info("Set thread count to: " + threadCount);
	}

	/**
	 * Clean Sage result folder by removing temporary files
	 */
	public static void cleanSageResultFolder(File folder, boolean keepResults) {
		if (!folder.exists()) {
			return;
		}
		
		File[] files = folder.listFiles();
		if (files != null) {
			for (File file : files) {
				if (file.isFile()) {
					// Keep important result files, remove temporary ones
					String name = file.getName();
					if (!keepResults || (!name.endsWith(".sage.tsv") && !name.endsWith(".sage.parquet"))) {
						if (name.endsWith(".bat") || name.endsWith(".json") || name.endsWith(".log")) {
							file.delete();
						}
					}
				}
			}
		}
	}

	public File getSageFile() {
		return sageFile;
	}

	public void run(int waitHours) {
		run(1, waitHours);
	}

	public void run(int threadPool, int waitHours) {

		if (this.parList.size() == 0) {
			return;
		}

		ExecutorService executor = Executors.newFixedThreadPool(threadPool);
		for (int i = 0; i < parList.size(); i++) {

			String pari = parList.get(i);

			LOGGER.info(pari);

			executor.submit(new Runnable() {

				@Override
				public void run() {
					// TODO Auto-generated method stub

					try {

						ArrayList<String> commands = new ArrayList<String>();
						commands.add("cmd.exe");
						commands.add("/c");
						String[] cs = pari.split("\t");
						for (String parString : cs) {
							commands.add(parString);
						}
						ProcessBuilder pb = new ProcessBuilder(commands);
						Process p = pb.start();
						BufferedInputStream in = new BufferedInputStream(p.getInputStream());
						BufferedReader inBr = new BufferedReader(new InputStreamReader(in));
						String line0 = "";
						String line1 = "";
						while ((line1 = inBr.readLine()) != null) {
							if (!line1.equals(line0)) {
								System.out.println(format.format(new Date()) + "\t" + line1);
								line0 = line1;
							}
						}

						if (p.waitFor() != 0) {
							if (p.exitValue() == 1)
								LOGGER.error("false");
						}
						inBr.close();
						in.close();

					} catch (IOException | InterruptedException e) {
						// TODO Auto-generated catch block
						LOGGER.error(e);
					}
				}
			});
		}

		try {

			executor.shutdown();

			boolean finish = executor.awaitTermination(waitHours, TimeUnit.HOURS);

			this.parList = new ArrayList<String>();

			if (finish) {
				System.out.println(format.format(new Date()) + "\t" + taskName + ": finished");
			} else {
				LOGGER.info(taskName
						+ ": task dosen't finish in a long time, please restart MetaLab after the Sage task finish");
				System.out.println(format.format(new Date()) + "\t" + taskName
						+ ": task dosen't finish in a long time, please restart MetaLab after the Sage task finish");
			}

		} catch (InterruptedException e) {
			// TODO Auto-generated catch block
			LOGGER.error(taskName + ": failed", e);
			System.out.println(format.format(new Date()) + "\t" + taskName + ": failed");
		}
	}

	public boolean runCMD(String cmd, JProgressBar bar) {

		try {

			int id = sageFile.getAbsolutePath().indexOf(":");
			if (id < 0) {
				if (id < 0) {
					LOGGER.error(taskName + ": cann't find the directory of sage.exe");
					System.out.println(
							format.format(new Date()) + "\t" + taskName + ": cann't find the directory of sage.exe");

					return false;
				}
			}

			PrintWriter writer = new PrintWriter(cmd);
			writer.println(sageFile.getAbsolutePath().substring(0, id) + ":");
			writer.println("cd " + sageFile.getParent());

			for (int i = 0; i < parList.size(); i++) {
				writer.println(sageFile.getName() + " " + parList.get(i));
			}

			writer.print("exit");
			writer.close();

		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			LOGGER.error(taskName + ": error in writing the executable command for sage.exe to " + cmd, e);
			System.out.println(format.format(new Date()) + "\t" + taskName
					+ ": error in writing the executable command for sage.exe to " + cmd);

			return false;
		}

		try {
			ArrayList<String> cmdList = new ArrayList<String>();
			cmdList.add("cmd.exe");
			cmdList.add("/c");
			cmdList.add(cmd);

			ProcessBuilder pb = new ProcessBuilder(cmdList);
			Process p = pb.start();

			BufferedInputStream in = new BufferedInputStream(p.getInputStream());
			BufferedReader inBr = new BufferedReader(new InputStreamReader(in));
			String lineStr;
			while ((lineStr = inBr.readLine()) != null) {
				LOGGER.info(lineStr);
				System.out.println(format.format(new Date()) + "\t" + taskName + ": " + lineStr);
			}

			if (p.waitFor() != 0) {
				if (p.exitValue() == 1)
					LOGGER.error("false");

				bar.setValue(100);
			}
			inBr.close();
			in.close();
		} catch (Exception e) {
			LOGGER.error(taskName + ": error in running the executable command by sage.exe", e);
			System.out.println(format.format(new Date()) + "\t" + taskName
					+ ": error in running the executable command by sage.exe");

			return false;
		}

		return true;
	}

}
