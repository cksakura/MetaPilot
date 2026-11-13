package uk.ac.quadram.metapilot.task.dda.mag.par;

import java.io.File;
import java.util.ArrayList;
import java.util.List;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.json.JSONArray;
import org.json.JSONObject;

import uk.ac.quadram.metapilot.core.msfragger.MSFraggerHeadlessRunner;
import uk.ac.quadram.metapilot.task.mag.GenomeSelectionCut;
import uk.ac.quadram.metapilot.core.opensearch.OpenSearchLayout;
import uk.ac.quadram.metapilot.task.par.BaseParameterIO;
import uk.ac.quadram.metapilot.task.par.MetaData;
import uk.ac.quadram.metapilot.task.par.SampleMetadata;

/**
 * Parameter persistence for the DDA MAG open-search workflow.
 * <p>
 * Uses the same sectioned JSON structure as the other MAG workflows —
 * {@code mag_database}, {@code input_output}, {@code enzyme_settings},
 * {@code additional_parameters}, {@code fixed_modifications},
 * {@code variable_modifications}, on top of {@code createJsonHeader(...)}'s
 * {@code workflow}/{@code version}/{@code generated_at}/{@code description}.
 * Workflow-specific settings go inside {@code additional_parameters}, exactly
 * as FragPipe's {@code search_mode} and {@code genome_selector_mode} do.
 * <p>
 * Experiments are read from and written to {@link MetaData} using
 * {@code file_path}/{@code experiment_name}/{@code fraction}. Missing
 * {@code fraction} defaults to 1 on load (same as pFind), so older GUI
 * exports without that field still produce a FlashLFQ-ready MetaData.
 * <p>
 * The filename comes from {@link OpenSearchLayout#PARAMETER_FILE} rather than a
 * duplicated literal, because {@code ResultFolderCleaner} detects the workflow
 * from that same name.
 */
public class DDAMagOpenSearchParameterIO extends BaseParameterIO {

	private static final Logger LOGGER = LogManager.getLogger(DDAMagOpenSearchParameterIO.class);

	private static final String LAST_PARAMETERS_FILENAME = "last_used_mag_opensearch_parameters.json";
	private static final String PREFERENCES_FILENAME = "mag_opensearch_preferences.json";

	@Override
	protected String getAnalysisParameterFileName() {
		return OpenSearchLayout.PARAMETER_FILE;
	}

	@Override
	protected String getLastUsedParametersFileName() {
		return LAST_PARAMETERS_FILENAME;
	}

	@Override
	protected String getPreferencesFileName() {
		return PREFERENCES_FILENAME;
	}

	@Override
	protected String getWorkflowType() {
		return "DDA_MAG_OPENSEARCH";
	}

	// ------------------------------------------------------------------ load

	public DDAMagOpenSearchParameters loadParameters(File inputFile) {
		try {
			LOGGER.info("Loading open-search parameters from: {}", inputFile.getAbsolutePath());
			JSONObject json = readJsonFromFile(inputFile);
			DDAMagOpenSearchParameters params = new DDAMagOpenSearchParameters();

			// --- input_output ------------------------------------------------
			if (json.has("input_output")) {
				JSONObject io = json.getJSONObject("input_output");
				params.setResult(io.optString("output_directory", ""));
				params.setMicroDb(io.optString("mag_database", ""));
				params.setMagDbCatalogueId(io.optString("mag_database", ""));
				params.setMagDbIdentifier(io.optString("mag_db_identifier", "MGYG"));
				params.setMagDbVersion(io.optString("mag_db_version", ""));
				params.setAppendHostDb(io.optBoolean("append_host_db", false));
				params.setHostDB(io.optString("host_database", ""));

				JSONArray experiments = io.optJSONArray("experiments");
				if (experiments != null) {
					String[] rawFiles = new String[experiments.length()];
					String[] expNames = new String[experiments.length()];
					int[] fractions = new int[experiments.length()];
					int[] replicates = new int[experiments.length()];
					List<String[]> rows = new ArrayList<>();
					for (int i = 0; i < experiments.length(); i++) {
						JSONObject exp = experiments.getJSONObject(i);
						String file = exp.optString("file_path", "");
						String name = exp.optString("experiment_name", "");
						int fraction = exp.optInt("fraction", 1);
						if (fraction <= 0) {
							fraction = 1;
						}
						rawFiles[i] = file;
						expNames[i] = name;
						fractions[i] = fraction;
						replicates[i] = 1;
						rows.add(new String[] { file, name, String.valueOf(fraction) });
					}
					// Same full MetaData reconstruction as DDAMagPfindParameterIO:
					// FlashLFQ ExperimentalDesign needs fractions/replicates; setting
					// only rawFiles/expNames leaves those arrays empty (length 0) and
					// mag-opensearch Step F crashes with ArrayIndexOutOfBoundsException.
					params.setMetadata(new MetaData(rawFiles, expNames, fractions, replicates));
					params.setExperimentData(rows);
				}
			}

			// --- mag_database ------------------------------------------------
			if (json.has("mag_database")) {
				JSONObject magDb = json.getJSONObject("mag_database");
				String catalogueId = magDb.optString("catalogue_id", "");
				if (!catalogueId.isBlank()) {
					params.setMagDbCatalogueId(catalogueId);
				}
				params.setMagDbIdentifier(magDb.optString("identifier", params.getMagDbIdentifier()));
				params.setMagDbVersion(magDb.optString("catalogue_version", params.getMagDbVersion()));
			}

			// --- enzyme_settings ---------------------------------------------
			if (json.has("enzyme_settings")) {
				JSONObject enzyme = json.getJSONObject("enzyme_settings");
				params.setEnzyme(migrateEnzymeName(enzyme,
						enzyme.optString("enzyme", MSFraggerHeadlessRunner.DEFAULT_ENZYME)));
				params.setMissCleavages(enzyme.optInt("missed_cleavages", 2));
				params.setDigestMode(enzyme.optString("digestion_mode", "Specific"));
			}

			// --- additional_parameters ---------------------------------------
			if (json.has("additional_parameters")) {
				JSONObject extra = json.getJSONObject("additional_parameters");
				params.setThreadCount(extra.optInt("threads", params.getThreadCount()));
				params.setMs2ScanMode(extra.optString("ms2_scan_mode", params.getMs2ScanMode()));
				params.setReduceDb(extra.optBoolean("generate_sample_specific_db", true));
				params.setSpectraCluster(extra.optBoolean("spectra_clustering", false));
				params.setGenomeSelectorMode(
						extra.optString("genome_selector_mode", params.getGenomeSelectorMode()));
				String storedStrategy = extra.optString("genome_selection_strategy", params.getGenomeSelectionStrategy());
				String strategyWarning = GenomeSelectionCut.strategyReplacementWarning(storedStrategy);
				if (strategyWarning != null) {
					LOGGER.warn(strategyWarning);
				}
				params.setGenomeSelectionStrategy(GenomeSelectionCut.strategyOrDefault(storedStrategy));
				GenomeSelectionCut.LoadedCoverage loaded = GenomeSelectionCut.readAdditionalKeys(extra,
						params.getGenomeSelectionPsmCoverage(), params.getGenomeSelectionStepGuard(),
						params.getGenomeSelectionMaxGenomes());
				params.setGenomeSelectionPsmCoverage(loaded.coverage);
				params.setGenomeSelectionStepGuard(loaded.stepGuard);
				params.setGenomeSelectionMaxGenomes(loaded.maxGenomes);
				if (loaded.replacementWarning != null) {
					LOGGER.warn(loaded.replacementWarning);
				}
				if (loaded.anyKey) {
					String warning = GenomeSelectionCut.unusedCoverageKeysWarning(params.getGenomeSelectionStrategy(),
							loaded.coverage, loaded.stepGuard, loaded.maxGenomes);
					if (warning != null) {
						LOGGER.warn(warning);
					}
				}
				params.setParallelTasks(extra.optInt("opensearch_parallel_tasks", params.getParallelTasks()));

				params.setSearchEngines(extra.optString("search_engines", params.getSearchEngines()));
				params.setMsfraggerMaxHeap(extra.optString("msfragger_max_heap", params.getMsfraggerMaxHeap()));
				params.setPrecursorMassLower(
						extra.optDouble("precursor_mass_lower", params.getPrecursorMassLower()));
				params.setPrecursorMassUpper(
						extra.optDouble("precursor_mass_upper", params.getPrecursorMassUpper()));
				params.setOutputReportTopN(extra.optInt("output_report_topn", params.getOutputReportTopN()));
				params.setRound1OutputReportTopN(
						extra.optInt("round1_output_report_topn", params.getRound1OutputReportTopN()));
				params.setRound1SearchMode(extra.optString("round1_search_mode", params.getRound1SearchMode()));
				params.setRound1ClosedPrecursorPpm(
						extra.optDouble("round1_closed_precursor_ppm", params.getRound1ClosedPrecursorPpm()));
				params.setClosedStageEnabled(extra.optBoolean("closed_stage_enabled", params.isClosedStageEnabled()));
				params.setClosedStagePrecursorPpm(
						extra.optDouble("closed_stage_precursor_ppm", params.getClosedStagePrecursorPpm()));
				params.setClosedOpenIntegration(
						extra.optString("closed_open_integration", params.getClosedOpenIntegration()));
				params.setClosedStageTopN(extra.optInt("closed_stage_topn", params.getClosedStageTopN()));
				params.setJointMs2Features(extra.optBoolean("joint_ms2_features", params.isJointMs2Features()));
				params.setCascadeMs2Features(extra.optBoolean("cascade_ms2_features", params.isCascadeMs2Features()));
				params.setCascadeClosedRerank(
						extra.optBoolean("cascade_closed_rerank", params.isCascadeClosedRerank()));
				params.setIdentificationVersion(
						extra.optString("identification_version", params.getIdentificationVersion()));
				// v2 changes two defaults; an explicit value in the file still wins.
				params.setSampleDecoyMethod(extra.optString("sample_decoy_method",
						params.isIdentificationV2() ? "pseudo_reverse" : params.getSampleDecoyMethod()));
				params.setJointMs2MaxOpenRank(extra.optInt("joint_ms2_max_open_rank",
						params.isIdentificationV2() ? Math.max(3, params.getJointMs2MaxOpenRank())
								: params.getJointMs2MaxOpenRank()));
				params.setJointFragmentTolerancePpm(
						extra.optDouble("joint_fragment_tolerance_ppm", params.getJointFragmentTolerancePpm()));
				params.setJointAcceptUnexplained(
						extra.optBoolean("joint_accept_unexplained", params.isJointAcceptUnexplained()));
				params.setJointCompareArms(extra.optBoolean("joint_compare_arms", params.isJointCompareArms()));
				params.setJointWithheldFeatures(
						extra.optString("joint_withheld_features", params.getJointWithheldFeatures()));
				params.setExpectedPeaksEnabled(
						extra.optBoolean("expected_peaks_enabled", params.isExpectedPeaksEnabled()));
				params.setExpectedPeakToleranceDa(
						extra.optDouble("expected_peak_tolerance_da", params.getExpectedPeakToleranceDa()));
				params.setExpectedPeakMinTargets(
						extra.optInt("expected_peak_min_targets", params.getExpectedPeakMinTargets()));
				params.setExpectedPeakMaxWindowFdr(
						extra.optDouble("expected_peak_max_window_fdr", params.getExpectedPeakMaxWindowFdr()));
				params.setExpectedPeakMinCorroboration(
						extra.optDouble("expected_peak_min_corroboration", params.getExpectedPeakMinCorroboration()));
				params.setExpectedPeakMaxAdmittedFraction(extra.optDouble("expected_peak_max_admitted_fraction",
						params.getExpectedPeakMaxAdmittedFraction()));
				params.setStage1RescoringEnabled(
						extra.optBoolean("stage1_rescoring_enabled", params.isStage1RescoringEnabled()));
				params.setStage1MsboosterEnabled(
						extra.optBoolean("stage1_msbooster_enabled", params.isStage1MsboosterEnabled()));
				params.setMsboosterMaxHeap(extra.optString("msbooster_max_heap", params.getMsboosterMaxHeap()));
				params.setMsboosterJarPath(extra.optString("msbooster_jar", params.getMsboosterJarPath()));
				params.setDiannExecutablePath(extra.optString("diann_executable", params.getDiannExecutablePath()));
				params.setUnimodOboPath(extra.optString("unimod_obo", params.getUnimodOboPath()));
				params.setDeltaMassEvidenceGateEnabled(
						extra.optBoolean("delta_mass_evidence_gate", params.isDeltaMassEvidenceGateEnabled()));
				params.setDeltaMassMaxNoGainFraction(
						extra.optDouble("delta_mass_max_no_gain_fraction", params.getDeltaMassMaxNoGainFraction()));
				params.setSubstitutionAnalysisEnabled(
						extra.optBoolean("substitution_analysis", params.isSubstitutionAnalysisEnabled()));
				params.setSubstitutionMassToleranceDa(
						extra.optDouble("substitution_mass_tolerance_da", params.getSubstitutionMassToleranceDa()));
				params.setSubstitutionMinFraction(
						extra.optDouble("substitution_min_fraction", params.getSubstitutionMinFraction()));
				params.setSubstitutionMinInDatabaseFraction(extra.optDouble(
						"substitution_min_in_db_fraction", params.getSubstitutionMinInDatabaseFraction()));
				params.setSequenceVariantCallingEnabled(extra.optBoolean("sequence_variant_calling",
						params.isSequenceVariantCallingEnabled()));
				params.setSequenceVariantMinSamples(
						extra.optInt("sequence_variant_min_samples", params.getSequenceVariantMinSamples()));
				params.getPangenomeVariantSettings().load(extra);
				params.getVariantFormRankingSettings().load(extra);
				params.setCrystalcEnabled(extra.optBoolean("crystalc_enabled", params.isCrystalcEnabled()));
				params.setCrystalcMaxHeap(extra.optString("crystalc_max_heap", params.getCrystalcMaxHeap()));
				params.setCrystalcJarPath(extra.optString("crystalc_jar", params.getCrystalcJarPath()));
				params.setCrystalcMinCharge(extra.optInt("crystalc_min_charge", params.getCrystalcMinCharge()));
				params.setCrystalcMaxCharge(extra.optInt("crystalc_max_charge", params.getCrystalcMaxCharge()));
				params.setCrystalcIsotopeNumber(
						extra.optInt("crystalc_isotope_number", params.getCrystalcIsotopeNumber()));
				params.setCrystalcPrecursorMassPpm(
						extra.optDouble("crystalc_precursor_mass_ppm", params.getCrystalcPrecursorMassPpm()));
				params.setCrystalcIsolationWindow(
						extra.optDouble("crystalc_isolation_window", params.getCrystalcIsolationWindow()));
				params.setZeroPeakSigmas(extra.optDouble("zero_peak_sigmas", params.getZeroPeakSigmas()));
				params.setZeroPeakSearchWindow(
						extra.optDouble("zero_peak_search_window", params.getZeroPeakSearchWindow()));
				params.setPeakMinR2(extra.optDouble("peak_min_r2", params.getPeakMinR2()));
				params.setMembershipFloorDa(extra.optDouble("membership_floor_da", params.getMembershipFloorDa()));
				params.setPeakMinArea(extra.optInt("peak_min_area", params.getPeakMinArea()));
				params.setMaxPeaksToFit(extra.optInt("max_peaks_to_fit", params.getMaxPeaksToFit()));
				params.setMaxFitIterations(extra.optInt("max_fit_iterations", params.getMaxFitIterations()));
				params.setPeakMinFitBins(extra.optInt("peak_min_fit_bins", params.getPeakMinFitBins()));
				params.setPeakMinSigmaBins(extra.optDouble("peak_min_sigma_bins", params.getPeakMinSigmaBins()));
				params.setPeakMinSnr(extra.optDouble("peak_min_snr", params.getPeakMinSnr()));
				params.setPeakMinSeparationBins(
						extra.optDouble("peak_min_separation_bins", params.getPeakMinSeparationBins()));
				params.setPeakMinSeparationSigmas(
						extra.optDouble("peak_min_separation_sigmas", params.getPeakMinSeparationSigmas()));
				params.setUseKdePeakDiscovery(
						extra.optBoolean("use_kde_peak_discovery", params.isUseKdePeakDiscovery()));
				params.setKdeBandwidthDa(extra.optDouble("kde_bandwidth_da", params.getKdeBandwidthDa()));
				params.setKdeScaleMultipliers(
						extra.optString("kde_scale_multipliers", params.getKdeScaleMultipliers()));
				params.setKdeMinPersistentScales(
						extra.optInt("kde_min_persistent_scales", params.getKdeMinPersistentScales()));
				params.setKdeMaxMembershipDa(
						extra.optDouble("kde_max_membership_da", params.getKdeMaxMembershipDa()));
				params.setKdeMinReproducibleFiles(
						extra.optInt("kde_min_reproducible_files", params.getKdeMinReproducibleFiles()));
				params.setKdeMinFileSupport(
						extra.optDouble("kde_min_file_support", params.getKdeMinFileSupport()));
				params.setKdeMaxModes(extra.optInt("kde_max_modes", params.getKdeMaxModes()));
				params.setKdeMaxSigmaMultiple(
						extra.optDouble("kde_max_sigma_multiple", params.getKdeMaxSigmaMultiple()));
				params.setKdeMinProminenceRatio(
						extra.optDouble("kde_min_prominence_ratio", params.getKdeMinProminenceRatio()));
				params.setRound2OutputReportTopN(
						extra.optInt("round2_output_report_topn", params.getRound2OutputReportTopN()));
				params.setRestrictPositivesToValidatedPeaks(extra.optBoolean(
						"restrict_positives_to_validated_peaks", params.isRestrictPositivesToValidatedPeaks()));
				params.setCombinedPsmInstances(
						extra.optBoolean("combined_psm_instances", params.isCombinedPsmInstances()));
				params.setMaxIsotopeOffset(extra.optInt("max_isotope_offset", params.getMaxIsotopeOffset()));
				params.setMinIsotopeFamilySatellites(extra.optInt("min_isotope_family_satellites",
						params.getMinIsotopeFamilySatellites()));
				params.setIsotopeMaxAbundanceIncrease(extra.optDouble("isotope_max_abundance_increase",
						params.getIsotopeMaxAbundanceIncrease()));
				params.setClassifierFolds(extra.optInt("classifier_folds", params.getClassifierFolds()));
				params.setClassifierMaxIterations(
						extra.optInt("classifier_max_iterations", params.getClassifierMaxIterations()));
				params.setClassifierTrainQValue(
						extra.optDouble("classifier_train_qvalue", params.getClassifierTrainQValue()));
				params.setLocalFdr(extra.optDouble("local_fdr", params.getLocalFdr()));
				params.setEntrapmentEnabled(extra.optBoolean("entrapment_enabled", params.isEntrapmentEnabled()));
				params.setEntrapmentSeed(extra.optLong("entrapment_seed", params.getEntrapmentSeed()));
				int entrapmentDigit = extra.optInt("entrapment_id_digit", params.getEntrapmentIdDigit());
				if (entrapmentDigit >= 1 && entrapmentDigit <= 9) {
					params.setEntrapmentIdDigit(entrapmentDigit);
				}
			}

			// --- modifications -----------------------------------------------
			if (json.has("fixed_modifications")) {
				JSONArray fixed = json.getJSONArray("fixed_modifications");
				String[] mods = new String[fixed.length()];
				for (int i = 0; i < fixed.length(); i++) {
					mods[i] = fixed.getString(i);
				}
				params.setFixMods(mods);
			}
			if (json.has("variable_modifications")) {
				JSONArray variable = json.getJSONArray("variable_modifications");
				String[] mods = new String[variable.length()];
				for (int i = 0; i < variable.length(); i++) {
					mods[i] = variable.getString(i);
				}
				params.setVariMods(mods);
			}

			List<String[]> metadataRows = SampleMetadata.fromJson(json);
			if (metadataRows != null) {
				params.setMetadataTableData(metadataRows);
				SampleMetadata.apply(params.getMetadata(), metadataRows);
			}
			return params;

		} catch (IllegalArgumentException e) {
			LOGGER.error("Rejected open-search parameters from {}: {}", inputFile.getAbsolutePath(), e.getMessage());
			throw e;
		} catch (Exception e) {
			LOGGER.error("Error loading open-search parameters from " + inputFile, e);
			return null;
		}
	}

	// ------------------------------------------------------------------ save

	public boolean saveParameters(DDAMagOpenSearchParameters params, File outputFile) {
		return saveParameters(params, outputFile, "Last used parameters");
	}

	public boolean saveParameters(DDAMagOpenSearchParameters params, File outputFile, String description) {
		try {
			File parentDir = outputFile.getParentFile();
			if (parentDir != null && !parentDir.exists() && !parentDir.mkdirs()) {
				LOGGER.error("Failed to create parent directory: {}", parentDir.getAbsolutePath());
				return false;
			}

			JSONObject json = createJsonHeader(description);

			// --- input_output ------------------------------------------------
			JSONObject ioSettings = new JSONObject();
			ioSettings.put("output_directory", params.getResult());
			ioSettings.put("mag_database",
					params.getMagDbCatalogueId() != null && !params.getMagDbCatalogueId().isBlank()
							? params.getMagDbCatalogueId()
							: params.getMicroDb());
			ioSettings.put("mag_db_identifier", params.getMagDbIdentifier());
			ioSettings.put("mag_db_version", params.getMagDbVersion());
			ioSettings.put("append_host_db", params.isAppendHostDb());
			ioSettings.put("host_database", params.getHostDB());

			JSONArray experiments = new JSONArray();
			MetaData metaData = params.getMetadata();
			if (metaData != null && metaData.getRawFiles() != null) {
				String[] rawFiles = metaData.getRawFiles();
				String[] expNames = metaData.getExpNames();
				int[] fractions = metaData.getFractions();
				for (int i = 0; i < rawFiles.length; i++) {
					JSONObject exp = new JSONObject();
					exp.put("file_path", rawFiles[i]);
					exp.put("experiment_name", expNames != null && i < expNames.length ? expNames[i] : "");
					int fraction = fractions != null && i < fractions.length && fractions[i] > 0 ? fractions[i] : 1;
					exp.put("fraction", fraction);
					experiments.put(exp);
				}
			} else if (params.getExperimentData() != null) {
				// Fall back to the panel-collected rows when metadata has not
				// been populated yet, so nothing is lost on an early save.
				for (String[] row : params.getExperimentData()) {
					JSONObject exp = new JSONObject();
					exp.put("file_path", row.length > 0 ? row[0] : "");
					exp.put("experiment_name", row.length > 1 ? row[1] : "");
					int fraction = 1;
					if (row.length > 2) {
						try {
							fraction = Integer.parseInt(row[2]);
							if (fraction <= 0) {
								fraction = 1;
							}
						} catch (NumberFormatException ignored) {
							fraction = 1;
						}
					}
					exp.put("fraction", fraction);
					experiments.put(exp);
				}
			}
			ioSettings.put("experiments", experiments);
			json.put("input_output", ioSettings);

			// --- mag_database ------------------------------------------------
			JSONObject magDbSettings = new JSONObject();
			String catalogueId = params.getMagDbCatalogueId();
			if ((catalogueId == null || catalogueId.isBlank()) && params.getMicroDb() != null) {
				catalogueId = params.getMicroDb();
			}
			magDbSettings.put("catalogue_id", catalogueId != null ? catalogueId : "");
			magDbSettings.put("identifier", params.getMagDbIdentifier());
			magDbSettings.put("catalogue_version",
					params.getMagDbVersion() != null ? params.getMagDbVersion() : "");
			json.put("mag_database", magDbSettings);

			// --- enzyme_settings ---------------------------------------------
			JSONObject enzymeSettings = new JSONObject();
			enzymeSettings.put("enzyme", params.getEnzyme());
			enzymeSettings.put(ENZYME_NAMING_KEY, ENZYME_NAMING_STANDARD);
			enzymeSettings.put("missed_cleavages", params.getMissCleavages());
			enzymeSettings.put("digestion_mode", params.getDigestModeString());
			json.put("enzyme_settings", enzymeSettings);

			// --- additional_parameters ---------------------------------------
			JSONObject additionalParams = new JSONObject();
			additionalParams.put("threads", params.getThreadCount());
			additionalParams.put("ms2_scan_mode", params.getMs2ScanMode());
			additionalParams.put("generate_sample_specific_db", params.isReduceDb());
			additionalParams.put("spectra_clustering", params.isSpectraCluster());
			additionalParams.put("genome_selector_mode", params.getGenomeSelectorMode());
			additionalParams.put("genome_selection_strategy", params.getGenomeSelectionStrategy());
			GenomeSelectionCut.writeAdditionalKeys(additionalParams, params.getGenomeSelectionPsmCoverage(),
					params.getGenomeSelectionStepGuard(), params.getGenomeSelectionMaxGenomes());
			additionalParams.put("opensearch_parallel_tasks", params.getParallelTasks());

			additionalParams.put("search_engines", params.getSearchEngines());
			additionalParams.put("msfragger_max_heap", params.getMsfraggerMaxHeap());
			additionalParams.put("precursor_mass_lower", params.getPrecursorMassLower());
			additionalParams.put("precursor_mass_upper", params.getPrecursorMassUpper());
			additionalParams.put("output_report_topn", params.getOutputReportTopN());
			additionalParams.put("round1_output_report_topn", params.getRound1OutputReportTopN());
			additionalParams.put("round1_search_mode", params.getRound1SearchMode());
			additionalParams.put("round1_closed_precursor_ppm", params.getRound1ClosedPrecursorPpm());
			additionalParams.put("closed_stage_enabled", params.isClosedStageEnabled());
			additionalParams.put("closed_stage_precursor_ppm", params.getClosedStagePrecursorPpm());
			additionalParams.put("closed_open_integration", params.getClosedOpenIntegration());
			additionalParams.put("closed_stage_topn", params.getClosedStageTopN());
			additionalParams.put("joint_ms2_features", params.isJointMs2Features());
			additionalParams.put("cascade_ms2_features", params.isCascadeMs2Features());
			additionalParams.put("cascade_closed_rerank", params.isCascadeClosedRerank());
			additionalParams.put("sample_decoy_method", params.getSampleDecoyMethod());
			additionalParams.put("identification_version", params.getIdentificationVersion());
			additionalParams.put("joint_ms2_max_open_rank", params.getJointMs2MaxOpenRank());
			additionalParams.put("joint_fragment_tolerance_ppm", params.getJointFragmentTolerancePpm());
			additionalParams.put("joint_accept_unexplained", params.isJointAcceptUnexplained());
			additionalParams.put("joint_compare_arms", params.isJointCompareArms());
			additionalParams.put("joint_withheld_features", params.getJointWithheldFeatures());
			additionalParams.put("expected_peaks_enabled", params.isExpectedPeaksEnabled());
			additionalParams.put("expected_peak_tolerance_da", params.getExpectedPeakToleranceDa());
			additionalParams.put("expected_peak_min_targets", params.getExpectedPeakMinTargets());
			additionalParams.put("expected_peak_max_window_fdr", params.getExpectedPeakMaxWindowFdr());
			additionalParams.put("expected_peak_min_corroboration", params.getExpectedPeakMinCorroboration());
			additionalParams.put("expected_peak_max_admitted_fraction", params.getExpectedPeakMaxAdmittedFraction());
			additionalParams.put("stage1_rescoring_enabled", params.isStage1RescoringEnabled());
			additionalParams.put("stage1_msbooster_enabled", params.isStage1MsboosterEnabled());
			additionalParams.put("msbooster_max_heap", params.getMsboosterMaxHeap());
			additionalParams.put("msbooster_jar", params.getMsboosterJarPath());
			additionalParams.put("diann_executable", params.getDiannExecutablePath());
			additionalParams.put("unimod_obo", params.getUnimodOboPath());
			additionalParams.put("delta_mass_evidence_gate", params.isDeltaMassEvidenceGateEnabled());
			additionalParams.put("delta_mass_max_no_gain_fraction", params.getDeltaMassMaxNoGainFraction());
			additionalParams.put("substitution_analysis", params.isSubstitutionAnalysisEnabled());
			additionalParams.put("substitution_mass_tolerance_da", params.getSubstitutionMassToleranceDa());
			additionalParams.put("substitution_min_fraction", params.getSubstitutionMinFraction());
			additionalParams.put("substitution_min_in_db_fraction", params.getSubstitutionMinInDatabaseFraction());
			additionalParams.put("sequence_variant_calling", params.isSequenceVariantCallingEnabled());
			additionalParams.put("sequence_variant_min_samples", params.getSequenceVariantMinSamples());
			params.getPangenomeVariantSettings().save(additionalParams);
			params.getVariantFormRankingSettings().save(additionalParams);
			additionalParams.put("crystalc_enabled", params.isCrystalcEnabled());
			additionalParams.put("crystalc_max_heap", params.getCrystalcMaxHeap());
			additionalParams.put("crystalc_jar", params.getCrystalcJarPath());
			additionalParams.put("crystalc_min_charge", params.getCrystalcMinCharge());
			additionalParams.put("crystalc_max_charge", params.getCrystalcMaxCharge());
			additionalParams.put("crystalc_isotope_number", params.getCrystalcIsotopeNumber());
			additionalParams.put("crystalc_precursor_mass_ppm", params.getCrystalcPrecursorMassPpm());
			additionalParams.put("crystalc_isolation_window", params.getCrystalcIsolationWindow());
			additionalParams.put("zero_peak_sigmas", params.getZeroPeakSigmas());
			additionalParams.put("zero_peak_search_window", params.getZeroPeakSearchWindow());
			additionalParams.put("peak_min_r2", params.getPeakMinR2());
			additionalParams.put("membership_floor_da", params.getMembershipFloorDa());
			additionalParams.put("peak_min_area", params.getPeakMinArea());
			additionalParams.put("max_peaks_to_fit", params.getMaxPeaksToFit());
			additionalParams.put("max_fit_iterations", params.getMaxFitIterations());
			additionalParams.put("peak_min_fit_bins", params.getPeakMinFitBins());
			additionalParams.put("peak_min_sigma_bins", params.getPeakMinSigmaBins());
			additionalParams.put("peak_min_snr", params.getPeakMinSnr());
			additionalParams.put("peak_min_separation_bins", params.getPeakMinSeparationBins());
			additionalParams.put("peak_min_separation_sigmas", params.getPeakMinSeparationSigmas());
			additionalParams.put("use_kde_peak_discovery", params.isUseKdePeakDiscovery());
			additionalParams.put("kde_bandwidth_da", params.getKdeBandwidthDa());
			additionalParams.put("kde_scale_multipliers", params.getKdeScaleMultipliers());
			additionalParams.put("kde_min_persistent_scales", params.getKdeMinPersistentScales());
			additionalParams.put("kde_max_membership_da", params.getKdeMaxMembershipDa());
			additionalParams.put("kde_min_reproducible_files", params.getKdeMinReproducibleFiles());
			additionalParams.put("kde_min_file_support", params.getKdeMinFileSupport());
			additionalParams.put("kde_max_modes", params.getKdeMaxModes());
			additionalParams.put("kde_max_sigma_multiple", params.getKdeMaxSigmaMultiple());
			additionalParams.put("kde_min_prominence_ratio", params.getKdeMinProminenceRatio());
			additionalParams.put("round2_output_report_topn", params.getRound2OutputReportTopN());
			additionalParams.put("restrict_positives_to_validated_peaks",
					params.isRestrictPositivesToValidatedPeaks());
			additionalParams.put("combined_psm_instances", params.isCombinedPsmInstances());
			additionalParams.put("max_isotope_offset", params.getMaxIsotopeOffset());
			additionalParams.put("min_isotope_family_satellites", params.getMinIsotopeFamilySatellites());
			additionalParams.put("isotope_max_abundance_increase", params.getIsotopeMaxAbundanceIncrease());
			additionalParams.put("classifier_folds", params.getClassifierFolds());
			additionalParams.put("classifier_max_iterations", params.getClassifierMaxIterations());
			additionalParams.put("classifier_train_qvalue", params.getClassifierTrainQValue());
			additionalParams.put("local_fdr", params.getLocalFdr());
			additionalParams.put("entrapment_enabled", params.isEntrapmentEnabled());
			additionalParams.put("entrapment_seed", params.getEntrapmentSeed());
			additionalParams.put("entrapment_id_digit", params.getEntrapmentIdDigit());
			json.put("additional_parameters", additionalParams);

			// --- modifications -----------------------------------------------
			JSONArray fixedMods = new JSONArray();
			if (params.getFixMods() != null) {
				for (String mod : params.getFixMods()) {
					fixedMods.put(mod);
				}
			}
			json.put("fixed_modifications", fixedMods);

			JSONArray variableMods = new JSONArray();
			if (params.getVariMods() != null) {
				for (String mod : params.getVariMods()) {
					variableMods.put(mod);
				}
			}
			json.put("variable_modifications", variableMods);
			json.put(SampleMetadata.JSON_KEY, SampleMetadata.toJson(params.getMetadataTableData()));

			boolean success = writeJsonToFile(json, outputFile);
			if (success) {
				LOGGER.info("Saved open-search parameters to: {}", outputFile.getAbsolutePath());
			}
			return success;

		} catch (Exception e) {
			LOGGER.error("Error saving open-search parameters to " + outputFile, e);
			return false;
		}
	}

	// --- BaseParameterIO contract -----------------------------------------

	@Override
	protected Object loadParametersFromFile(File inputFile) {
		return loadParameters(inputFile);
	}

	@Override
	protected boolean saveParametersToFile(Object parameters, File outputFile, String description) {
		if (!(parameters instanceof DDAMagOpenSearchParameters)) {
			LOGGER.error("Invalid parameter type: expected DDAMagOpenSearchParameters, got {}",
					parameters == null ? "null" : parameters.getClass().getName());
			return false;
		}
		return saveParameters((DDAMagOpenSearchParameters) parameters, outputFile, description);
	}

	@Override
	protected Object createDefaultParameters() {
		return new DDAMagOpenSearchParameters();
	}
	/**
	 * Marks a parameter file as using the corrected enzyme naming.
	 * <p>
	 * Before the correction, {@code "Trypsin"} in these files meant MSFragger's
	 * {@code stricttrypsin} - cleavage after K/R including before proline - and
	 * {@code "Trypsin/P"} meant the conventional proline-restricted digest. The
	 * names now mean what they say. A file written before the change therefore
	 * asks for a different digest than its author intended, and reloading one
	 * would silently re-run a study with a different enzyme.
	 * <p>
	 * The shared header's {@code version} field is hard-coded and never read
	 * back, so it cannot date a file. The presence of this key can: it is written
	 * from now on and absent from everything older.
	 */
	private static final String ENZYME_NAMING_KEY = "enzyme_naming";
	private static final String ENZYME_NAMING_STANDARD = "standard";

	/**
	 * Translates a pre-correction enzyme name to the one that means the same
	 * digest today. Only the two inverted pairs move; everything else is
	 * unchanged, and a file already carrying {@link #ENZYME_NAMING_KEY} is left
	 * exactly as written.
	 */
	private static String migrateEnzymeName(JSONObject enzymeSettings, String enzyme) {
		if (enzymeSettings.has(ENZYME_NAMING_KEY) || enzyme == null) {
			return enzyme;
		}
		String migrated;
		switch (enzyme.trim().toLowerCase(java.util.Locale.ROOT)) {
		case "trypsin":
			migrated = "Trypsin/P";
			break;
		case "trypsin/p":
		case "trypsin_p":
			migrated = "Trypsin";
			break;
		case "lysc":
			migrated = "LysC/P";
			break;
		case "lysc/p":
		case "lysc_p":
			migrated = "LysC";
			break;
		default:
			return enzyme;
		}
		LOGGER.warn("This parameter file predates the enzyme-naming correction, where \"{}\" meant the digest now "
				+ "called \"{}\". Loading it as \"{}\" so the study runs the enzyme it was written with. Re-save the "
				+ "file to make this permanent.", enzyme, migrated, migrated);
		return migrated;
	}

}
