package uk.ac.quadram.metapilot.core.runner;

import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Locale;
import java.util.Map;

import uk.ac.quadram.metapilot.core.mag.DdaMagPostAnalysisService;
import uk.ac.quadram.metapilot.core.mag.GenomeEvidenceTest;
import uk.ac.quadram.metapilot.core.opensearch.OpenSearchMetaPeptide;
import uk.ac.quadram.metapilot.core.opensearch.OpenSearchMetaPsm;
import uk.ac.quadram.metapilot.core.opensearch.OpenSearchPsm;
import uk.ac.quadram.metapilot.core.prodb.TargetDecoyFasta;
import uk.ac.quadram.metapilot.core.workflow.ProgressCallback;
import uk.ac.quadram.metapilot.quant.flashLFQ.FlashLfqQuanPepReader;
import uk.ac.quadram.metapilot.quant.flashLFQ.FlashLfqQuanPeptide;
import uk.ac.quadram.metapilot.task.io.pep.MetaPSM;
import uk.ac.quadram.metapilot.task.io.pep.MetaPeptide;
import uk.ac.quadram.metapilot.task.io.qpx.QpxPsmExport;
import uk.ac.quadram.metapilot.task.io.pro.MetaProtein;
import uk.ac.quadram.metapilot.task.mag.MagDbItem;
import uk.ac.quadram.metapilot.task.par.MetaConstants;

/**
 * Post-analysis export for the open-search workflow, mirroring
 * {@link FragpipePostAnalysisTsvExporter}.
 * <p>
 * This is the seam between the engine-specific part of this workflow and the
 * shared MAG skeleton: it aggregates validated PSMs into {@link MetaPeptide}s,
 * hands them to {@link DdaMagPostAnalysisService} — the same generic service
 * every DDA MAG engine uses for genome scoring, protein inference and razor
 * assignment — and writes the standard deliverables through
 * {@link DdaMagReportExportService}.
 * <p>
 * Two differences from the FragPipe exporter, both consequences of this
 * workflow's design rather than choices:
 * <ul>
 * <li>The decoy symbol is {@code rev_}, not {@code shuf_}. FragPipe needs
 * {@code shuf_} because Philosopher strips {@code rev_} from its quantified
 * peptide list; this workflow emits its own lists with decoys retained.</li>
 * <li>Peptides are keyed on the <em>modified</em> sequence including the delta
 * mass, so the same backbone bearing different modifications aggregates as
 * distinct peptides — which is the point of an open search.</li>
 * </ul>
 * <h2>Identification vs. quantification</h2>
 * {@link #aggregatePeptides} builds the identification truth for every
 * peptide — score, decoy flag, delta mass/peak id, localized sequence,
 * charges, MS2 counts, and a first protein set from the validated PSMs alone
 * — with zero intensity, because a PSM list has no quantification in it.
 * When a FlashLFQ {@code QuantifiedPeptides.tsv} is supplied (see
 * {@code quantifiedPeptidesTsv} below), {@link #rematchWithFlashLfq} then
 * overlays FlashLFQ's real per-sample intensities onto that same peptide set
 * and unions in FlashLFQ's own protein-group assignment, which was resolved
 * from every sample's evidence in one combined run rather than one sample's
 * own psms.tsv rows. Without it, every peptide keeps phase-1's zero
 * intensity and its psms-only protein set, exactly as before this existed.
 */
public final class OpenSearchPostAnalysisTsvExporter {

	private OpenSearchPostAnalysisTsvExporter() {
	}

	public static final class Result {
		public File finalPeptidesTsv;
		public File razorProteinsTsv;
		public File allProteinsTsv;
		public int peptideCount;
		public int proteinCount;
		public HashMap<String, Double> genomeScoreMap;

		// The fields below mirror FragpipePostAnalysisTsvExporter.Result exactly -
		// same names, same source (DdaMagPostAnalysisService.Result, the shared
		// service both exporters call) - so MagTaxonomyAnalysisService and
		// MagFunctionalAnnotationService can be wired up here the same way
		// MagFragPipeHeadlessRunner already wires them up from that exporter's
		// Result, rather than each engine inventing its own Request-building code.
		public String[] intensityTitles;
		public MetaProtein[] metaProteins;
		public HashMap<String, double[]> genomeIntensityMap;
		/** genome -> double[]{S, p_emp, q_value} when PSM-score/top-k scoring is used; null otherwise */
		public HashMap<String, double[]> genomeDiaStats;
		/** The genome test behind {@link #genomeDiaStats} when it was {@link GenomeEvidenceTest}; else null. */
		public GenomeEvidenceTest.Result genomeEvidence;
		public HashMap<String, Integer> psmCountMap;
		public HashMap<String, HashSet<String>> totalPepCountMap;
		public HashMap<String, HashSet<String>> razorPepCountMap;
		public HashMap<String, HashSet<String>> uniquePepCountMap;
		public HashMap<String, HashSet<String>> totalPepWithILMap;
		public HashMap<String, HashSet<String>> razorPepWithILMap;
		public HashMap<String, HashSet<String>> uniquePepWithILMap;
		/** See {@link DdaMagPostAnalysisService.Result#peptidesWithoutConfidentRazorProtein}. */
		public int peptidesWithoutConfidentRazorProtein;
	}

	/**
	 * @param magDb        selected MAG database
	 * @param validated    validated PSMs from {@code ModificationValidator}
	 * @param sampleNames  sample base names, in the order used for intensity columns
	 * @param outDir       where the deliverables land, and the search root
	 *                     whose stage-2 parameters give the fixed modifications
	 */
	public static Result exportFinalTsvs(MagDbItem magDb, List<OpenSearchPsm> validated, List<String> sampleNames,
			File outDir, File spectraDir, String taskName, ProgressCallback cb) throws Exception {
		return exportFinalTsvs(magDb, validated, sampleNames, outDir, outDir, spectraDir, taskName, null, cb);
	}

	/**
	 * @param outDir                where the deliverables land: an engine
	 *                              configuration's {@code OpenSearchLayout.modeRoot},
	 *                              with {@code id.tsv} in its {@code combined/}
	 * @param searchRoot            the engine search root holding
	 *                              {@code search_results/}, read for the stage-2
	 *                              fixed modifications
	 * @param quantifiedPeptidesTsv FlashLFQ's own {@code QuantifiedPeptides.tsv}
	 *                              from quantifying {@code validated_psms.tsv},
	 *                              or null to keep phase-1 behaviour (identify
	 *                              only, zero intensity). See
	 *                              {@link #rematchWithFlashLfq} for what
	 *                              supplying it changes.
	 */
	public static Result exportFinalTsvs(MagDbItem magDb, List<OpenSearchPsm> validated, List<String> sampleNames,
			File outDir, File searchRoot, File spectraDir, String taskName, File quantifiedPeptidesTsv,
			ProgressCallback cb) throws Exception {
		return exportFinalTsvs(magDb, validated, sampleNames, outDir, searchRoot, spectraDir, taskName,
				quantifiedPeptidesTsv, null, cb);
	}

	/**
	 * @param genomeEvidence genome and protein statistics from
	 *                       {@link GenomeEvidenceTest}, or null for the
	 *                       resampling genome test on the validated decoys
	 */
	public static Result exportFinalTsvs(MagDbItem magDb, List<OpenSearchPsm> validated, List<String> sampleNames,
			File outDir, File searchRoot, File spectraDir, String taskName, File quantifiedPeptidesTsv,
			GenomeEvidenceTest.Result genomeEvidence, ProgressCallback cb) throws Exception {
		if (cb == null) {
			cb = new ProgressCallback.NoOp();
		}
		if (magDb == null) {
			throw new IllegalArgumentException("magDb is null");
		}
		if (outDir == null) {
			throw new IllegalArgumentException("outDir is null");
		}
		outDir.mkdirs();

		HashMap<String, Integer> fileIdMap = new HashMap<>();
		ArrayList<String> fileList = new ArrayList<>();
		for (String sample : sampleNames) {
			if (sample == null || sample.isBlank() || fileIdMap.containsKey(sample)) {
				continue;
			}
			fileIdMap.put(sample, fileIdMap.size());
			fileList.add(sample);
		}
		String[] intensityTitles = fileList.toArray(new String[0]);

		// Decoys remain in validated_psms.tsv/scored_psms.tsv because they are the
		// audit trail for FDR, but they are not biological identifications and must
		// never enter id.tsv or any peptide/protein/taxon report. Their one other
		// use is as the genome-scoring null (Input.genomeBackgroundPsms): the
		// target list carries no decoy-only peptides, so without them that null
		// would be empty.
		List<OpenSearchPsm> targetPsms = targetsOnly(validated);
		if (targetPsms.isEmpty()) {
			throw new IllegalArgumentException("No validated target PSMs are available for post-analysis");
		}
		List<OpenSearchPsm> decoyPsms = genomeEvidence == null ? decoysOnly(validated) : new ArrayList<>();
		cb.onMessage(taskName + ": aggregating " + targetPsms.size() + " validated target PSMs into peptides"
				+ (genomeEvidence == null
						? " (" + decoyPsms.size() + " decoys used only as the genome-scoring background) ..."
						: " (genomes already tested on the spectrum winners) ..."));
		MetaPeptide[] peptides = aggregatePeptides(targetPsms, fileIdMap);
		MetaPeptide[] decoyPeptides = aggregatePeptides(decoyPsms, fileIdMap);
		if (quantifiedPeptidesTsv != null && quantifiedPeptidesTsv.isFile() && quantifiedPeptidesTsv.length() > 0) {
			Map<String, FlashLfqQuanPeptide> quantified = readQuantifiedPeptides(quantifiedPeptidesTsv, fileList);
			rematchWithFlashLfq(peptides, quantified, taskName + ": FlashLFQ rematch", cb);
			// Genome scoring weights each peptide by its intensity; an unquantified
			// decoy would sit at the intensity floor and weaken the null.
			if (decoyPeptides.length > 0) {
				rematchWithFlashLfq(decoyPeptides, quantified,
						taskName + ": FlashLFQ rematch of genome-scoring decoys", cb);
			}
		} else {
			cb.onMessage(taskName + ": no FlashLFQ quantification supplied; peptides keep phase-1 zero intensity");
		}
		Map<String, Double> staticMods = searchStaticMods(searchRoot, fileList, cb);
		MetaPSM[] psms = toMetaPsms(targetPsms, staticMods);
		MetaPSM[] decoyMetaPsms = toMetaPsms(decoyPsms, staticMods);
		cb.onMessage(taskName + ": peptides=" + peptides.length + ", psms=" + psms.length);

		DdaMagPostAnalysisService.Result analysis = DdaMagPostAnalysisService.analyze(
				new DdaMagPostAnalysisService.Input(peptides, psms, fileIdMap, fileList,
						TargetDecoyFasta.DECOY_TAG, magDb.getIdentifier(),
						OpenSearchPostAnalysisTsvExporter::resolveFileIndex, true, decoyPeptides, decoyMetaPsms,
						genomeEvidence));
		String[] proteinScoreColumns = genomeEvidence != null ? DdaMagReportExportService.PICKED_PROTEIN_SCORE_COLUMNS
				: DdaMagReportExportService.RESAMPLING_PROTEIN_SCORE_COLUMNS;

		cb.onMessage(taskName + ": genome scoring produced " + analysis.genomeScoreMap.size() + " genomes");

		File combinedDir = uk.ac.quadram.metapilot.core.opensearch.OpenSearchLayout.combinedDir(outDir);
		combinedDir.mkdirs();
		File pepOut = new File(outDir, "final_peptides.tsv");
		File proOut = new File(outDir, "razor_proteins.tsv");
		File allProOut = new File(outDir, "all_proteins.tsv");

		cb.onMessage(taskName + ": writing id.tsv");
		DdaMagReportExportService.writeIdentificationSummary(combinedDir, outDir, spectraDir, fileIdMap,
				analysis.filePsmCount, analysis.fileSequenceSets, null);

		cb.onMessage(taskName + ": writing " + pepOut.getName());
		DdaMagReportExportService.writePeptideReport(pepOut, analysis.peptideList, analysis.psmMap,
				analysis.pepPsmCountMap, analysis.pepChargeMap, analysis.proteinIdMap, fileList,
				MetaConstants.labelFree, intensityTitles, buildRazorProteinMap(peptides, analysis.razorProteins));

		cb.onMessage(taskName + ": writing " + proOut.getName());
		DdaMagReportExportService.writeProteinReport(proOut, analysis.proteins, analysis.proteinIntensityMap,
				analysis.proteinTotalIntenMap, analysis.totalPepCountMap, analysis.razorPepCountMap,
				analysis.psmCountMap, MetaConstants.labelFree, intensityTitles, analysis.peptideList,
				analysis.allProteinScoreMap, analysis.psmMap, proteinScoreColumns);

		cb.onMessage(taskName + ": writing " + allProOut.getName());
		DdaMagReportExportService.writeAllProteinReport(allProOut, analysis.allProteins,
				analysis.allProteinIntensityMap, analysis.allProteinTotalIntenMap, analysis.allProteinPepCountMap,
				analysis.razorPepCountMap, analysis.allProteinPsmCountMap, MetaConstants.labelFree, intensityTitles,
				analysis.peptideList, analysis.allProteinScoreMap, analysis.psmMap, proteinScoreColumns);

		writeQpxPsms(validated, staticMods, outDir, cb);

		Result result = new Result();
		result.finalPeptidesTsv = pepOut;
		result.razorProteinsTsv = proOut;
		result.allProteinsTsv = allProOut;
		result.peptideCount = peptides.length;
		result.proteinCount = analysis.proteins == null ? 0 : analysis.proteins.length;
		result.genomeScoreMap = analysis.genomeScoreMap;
		result.intensityTitles = intensityTitles;
		result.metaProteins = analysis.metaProteins;
		result.genomeIntensityMap = analysis.genomeIntensityMap;
		result.genomeDiaStats = analysis.genomeDiaStats;
		result.genomeEvidence = genomeEvidence;
		result.psmCountMap = analysis.psmCountMap;
		result.totalPepCountMap = analysis.totalPepCountMap;
		result.razorPepCountMap = analysis.razorPepCountMap;
		result.uniquePepCountMap = analysis.uniquePepCountMap;
		result.totalPepWithILMap = analysis.totalPepWithILMap;
		result.razorPepWithILMap = analysis.razorPepWithILMap;
		result.uniquePepWithILMap = analysis.uniquePepWithILMap;
		result.peptidesWithoutConfidentRazorProtein = analysis.peptidesWithoutConfidentRazorProtein;
		return result;
	}

	/**
	 * Writes the PSM layer as QPX psm.parquet beside final_peptides.tsv, on its
	 * own, for a run that keeps its final tables. Never throws.
	 *
	 * @param searchRoot the engine search root, read for the stage-2 fixed
	 *                   modifications
	 * @param outDir     the folder holding final_peptides.tsv
	 */
	public static File writeQpxPsms(List<OpenSearchPsm> validated, List<String> sampleNames, File searchRoot,
			File outDir, ProgressCallback cb) {
		return writeQpxPsms(validated, searchStaticMods(searchRoot, sampleNames, cb), outDir, cb);
	}

	/**
	 * Unlike the peptide layer this keeps the validated decoys: they are the
	 * population the PSM q-values were estimated against, which QPX carries via
	 * is_decoy.
	 */
	private static File writeQpxPsms(List<OpenSearchPsm> validated, Map<String, Double> staticMods, File outDir,
			ProgressCallback cb) {
		List<MetaPSM> qpxPsms = new ArrayList<>(validated.size());
		for (OpenSearchPsm psm : validated) {
			qpxPsms.add(new OpenSearchMetaPsm(psm, staticMods));
		}
		return QpxPsmExport.writeQuietly(outDir, qpxPsms, cb);
	}

	/**
	 * The fixed modifications of the stage-2 search, from the first sample's
	 * fragger.params (the same file SageOpenSearchRunner mirrors into Sage, so both
	 * engines share them). Empty - and the QPX peptidoform then omits them - if no
	 * params file is found.
	 */
	static Map<String, Double> searchStaticMods(File resultRoot, List<String> samples, ProgressCallback cb) {
		if (resultRoot != null && samples != null) {
			for (String sample : samples) {
				File params = uk.ac.quadram.metapilot.core.opensearch.OpenSearchLayout.sampleDbParams(resultRoot,
						sample);
				if (!params.isFile()) {
					continue;
				}
				try {
					uk.ac.quadram.metapilot.core.sage.SageOpenSearchRunner.Settings settings =
							new uk.ac.quadram.metapilot.core.sage.SageOpenSearchRunner.Settings();
					uk.ac.quadram.metapilot.core.sage.SageOpenSearchRunner.applyMsfraggerParams(settings, params);
					return new LinkedHashMap<>(settings.staticMods);
				} catch (Exception e) {
					cb.onMessage("QPX PSM export: could not read fixed modifications from " + params + ": "
							+ e.getMessage());
				}
			}
		}
		cb.onMessage("QPX PSM export: no stage-2 fragger.params found; fixed modifications omitted from peptidoforms");
		return new LinkedHashMap<>();
	}

	static List<OpenSearchPsm> targetsOnly(List<OpenSearchPsm> psms) {
		List<OpenSearchPsm> targets = new ArrayList<>();
		if (psms != null) {
			for (OpenSearchPsm psm : psms) {
				if (psm != null && psm.isTarget()) {
					targets.add(psm);
				}
			}
		}
		return targets;
	}

	static List<OpenSearchPsm> decoysOnly(List<OpenSearchPsm> psms) {
		List<OpenSearchPsm> decoys = new ArrayList<>();
		if (psms != null) {
			for (OpenSearchPsm psm : psms) {
				if (psm != null && !psm.isTarget()) {
					decoys.add(psm);
				}
			}
		}
		return decoys;
	}

	private static MetaPSM[] toMetaPsms(List<OpenSearchPsm> psms, Map<String, Double> staticMods) {
		MetaPSM[] out = new MetaPSM[psms.size()];
		for (int i = 0; i < out.length; i++) {
			out[i] = new OpenSearchMetaPsm(psms.get(i), staticMods);
		}
		return out;
	}

	/**
	 * Groups PSMs by modified sequence. The delta mass is part of the key, so a
	 * modified and unmodified form of the same backbone stay separate — the
	 * distinction this workflow exists to make.
	 */
	private static MetaPeptide[] aggregatePeptides(List<OpenSearchPsm> validated,
			HashMap<String, Integer> fileIdMap) {
		Map<String, List<OpenSearchPsm>> byModSeq = new LinkedHashMap<>();
		for (OpenSearchPsm psm : validated) {
			byModSeq.computeIfAbsent(canonicalModifiedSequenceOf(psm), k -> new ArrayList<>()).add(psm);
		}

		int fileCount = Math.max(1, fileIdMap.size());
		List<MetaPeptide> peptides = new ArrayList<>(byModSeq.size());

		for (Map.Entry<String, List<OpenSearchPsm>> entry : byModSeq.entrySet()) {
			List<OpenSearchPsm> group = entry.getValue();
			OpenSearchPsm best = null;
			for (OpenSearchPsm psm : group) {
				if (best == null || score(psm) > score(best)) {
					best = psm;
				}
			}

			HashSet<String> proteinSet = new HashSet<>();
			HashSet<Integer> chargeSet = new HashSet<>();
			int[] ms2Counts = new int[fileCount];
			for (OpenSearchPsm psm : group) {
				// The FULL protein list, not the delegate. Razor assignment and
				// genome scoring downstream are entirely about peptides that map
				// to more than one protein; with one protein per peptide every
				// peptide looks trivially unique and razor inference degenerates.
				// Matches FragpipePSMTsvReader, which merges "Protein" with
				// "Mapped Proteins".
				for (String protein : psm.getAllProteins()) {
					if (protein != null && !protein.isBlank()) {
						proteinSet.add(protein.trim());
					}
				}
				if (proteinSet.isEmpty() && psm.getProtein() != null && !psm.getProtein().isBlank()) {
					proteinSet.add(psm.getProtein().trim());
				}
				chargeSet.add(psm.getCharge());
				Integer fileId = fileIdMap.get(psm.getFileName());
				if (fileId != null && fileId < ms2Counts.length) {
					ms2Counts[fileId]++;
				}
			}

			// Phase 1 has no quantification; the report writers still expect an
			// intensity array of the right width.
			double[] intensity = new double[fileCount];

			OpenSearchMetaPeptide peptide = new OpenSearchMetaPeptide(best.getPeptide(), entry.getKey(),
					proteinSet.toArray(new String[0]), intensity, new int[fileCount],
					best.isClosedStage() ? best.getClosedModificationMass() : best.getNormalizedMassDiff(),
					best.getPeakId(), best.getLocalizedSequence());

			peptide.setScore(score(best));
			peptide.setDecoy(!best.isTarget());
			peptide.setMs2Counts(ms2Counts);
			peptide.setTotalMS2Count(group.size());

			int[] charges = new int[chargeSet.size()];
			int i = 0;
			for (Integer charge : chargeSet) {
				charges[i++] = charge;
			}
			peptide.setCharges(charges);

			peptides.add(peptide);
		}
		return peptides.toArray(new MetaPeptide[0]);
	}

	/**
	 * Delegates to {@link OpenSearchPsm#getCanonicalModifiedSequence()}, the
	 * single source of truth for this string. Kept as a static here because
	 * every existing caller (this class, {@code MagOpenSearchHeadlessRunner})
	 * already calls it this way; moving the logic to {@code OpenSearchPsm}
	 * itself was what let {@code OpenSearchMetaPsm} - a different package -
	 * reuse it too, so the PSM-side and peptide-side "modSeq" join key are
	 * always the same string. See that method's javadoc for why the two used
	 * to diverge.
	 */
	static String canonicalModifiedSequenceOf(OpenSearchPsm psm) {
		return psm.getCanonicalModifiedSequence();
	}

	/**
	 * Resolves {@link OpenSearchPsm#getCanonicalModifiedSequence()} for every
	 * supplied PSM, exactly once, before anything reads it. Must run first -
	 * {@code MagOpenSearchHeadlessRunner} calls this immediately after
	 * {@code ModificationValidator} returns, before {@code writeValidatedPsms},
	 * and before {@link #aggregatePeptides} or {@code OpenSearchMetaPsm} ever
	 * see these PSMs - so the "Full Sequence" column FlashLFQ reads, the
	 * peptide grouping key {@link #aggregatePeptides} builds, and
	 * {@code OpenSearchMetaPsm}'s own join key are all the same string,
	 * decided the same way.
	 * <p>
	 * Groups PSMs by {@link OpenSearchPsm#getModificationFamilyKey()} -
	 * peptide backbone plus delta-mass family, never a specific localized
	 * residue - then, within each group, picks the single highest-scoring PSM
	 * using the same tie-break {@link #score} and {@link #aggregatePeptides}
	 * already use to pick their own representative. That one PSM's own
	 * {@link OpenSearchPsm#getLocalizedSequence()} decides where the
	 * {@code [+mass]} bracket goes - immediately after the first residue it
	 * marks lowercase, or at the end of the sequence when even the
	 * representative has no localized site - and every member of the group,
	 * not just the representative, is stamped with that one resulting string.
	 * <p>
	 * This is deliberate, not an approximation: {@link OpenSearchPsm#getCanonicalModifiedSequence()}'s
	 * javadoc explains why a peptide's identity must never depend on which
	 * spectrum's localization guess happened to produce it. Once grouped this
	 * way, which single site gets displayed is cosmetic - it does not change
	 * which rows this peptide's PSMs, or its FlashLFQ intensity, are split
	 * across.
	 */
	static void resolveCanonicalModifiedSequences(List<OpenSearchPsm> psms) {
		if (psms == null || psms.isEmpty()) {
			return;
		}
		Map<String, List<OpenSearchPsm>> byFamily = new LinkedHashMap<>();
		for (OpenSearchPsm psm : psms) {
			byFamily.computeIfAbsent(psm.getModificationFamilyKey(), k -> new ArrayList<>()).add(psm);
		}
		for (List<OpenSearchPsm> group : byFamily.values()) {
			OpenSearchPsm best = null;
			for (OpenSearchPsm psm : group) {
				if (best == null || score(psm) > score(best)) {
					best = psm;
				}
			}
			String resolvedSequence = decorateWithRepresentativeSite(best);
			double resolvedMass = representativeMonoisotopicMass(best);
			int resolvedSite = representativeSiteIndex(best);
			for (OpenSearchPsm psm : group) {
				psm.setResolvedModifiedSequence(resolvedSequence);
				psm.setResolvedMonoisotopicMass(resolvedMass);
				psm.setResolvedSiteIndex(resolvedSite);
			}
		}
	}

	/**
	 * The delta mass both {@link #decorateWithRepresentativeSite} and
	 * {@link #representativeMonoisotopicMass} key off, computed once so they
	 * can never disagree: the fitted peak's centroid, falling back to the
	 * representative's own normalized mass diff when no centroid was fitted.
	 * {@link Double#NaN} when the representative has no peak at all
	 * ({@code peakId < 0}) - the isotope-zero collapse (|mass| &lt;= 0.05) is
	 * a range check each caller applies itself, not encoded here, because the
	 * two callers do different things once they know it applies.
	 */
	private static double canonicalMassOf(OpenSearchPsm representative) {
		if (representative.getPeakId() < 0) {
			return Double.NaN;
		}
		double canonicalMass = representative.getPeakCentroid();
		if (Double.isNaN(canonicalMass)) {
			canonicalMass = representative.getNormalizedMassDiff();
		}
		return canonicalMass;
	}

	/**
	 * Builds one group's final {@code getCanonicalModifiedSequence()} string
	 * from its representative PSM: the family key
	 * ({@code representative.getModificationFamilyKey()}) with the
	 * {@code [+mass]} bracket relocated to sit immediately after the
	 * representative's own localized residue, when it has one. Peptide
	 * letters stay upper case throughout - only the bracket's position
	 * changes, never its numeric content or the surrounding letters' case.
	 */
	/**
	 * The 0-based residue {@link #decorateWithRepresentativeSite} puts the bracket on,
	 * or -1 when it returns the position-free family key. Same conditions, kept next
	 * to it so the two cannot drift.
	 */
	static int representativeSiteIndex(OpenSearchPsm representative) {
		if (representative.isClosedStage()) {
			double[] mods = representative.getClosedModifications();
			for (int i = 1; mods != null && i < mods.length; i++) {
				if (mods[i] != 0.0) {
					return i - 1;
				}
			}
			return -1;
		}
		double canonicalMass = canonicalMassOf(representative);
		if (Double.isNaN(canonicalMass) || Math.abs(canonicalMass) <= 0.05) {
			return -1;
		}
		String localized = representative.getLocalizedSequence();
		String peptide = representative.getPeptide();
		int site = representative.getFirstLocalizedSiteIndex();
		if (site < 0 || localized == null || peptide == null || localized.length() != peptide.length()) {
			return -1;
		}
		return Math.min(site, peptide.length() - 1);
	}

	private static String decorateWithRepresentativeSite(OpenSearchPsm representative) {
		String familyKey = representative.getModificationFamilyKey();
		if (representative.isClosedStage()) {
			return familyKey;
		}
		// The family key already carries no position for the unexplained and
		// isotope-zero branches (see getModificationFamilyKey()); only a real,
		// in-peak modification has a residue worth showing.
		double canonicalMass = canonicalMassOf(representative);
		if (Double.isNaN(canonicalMass) || Math.abs(canonicalMass) <= 0.05) {
			return familyKey;
		}
		String localized = representative.getLocalizedSequence();
		String peptide = representative.getPeptide();
		int site = representative.getFirstLocalizedSiteIndex();
		if (site < 0 || localized == null || peptide == null || localized.length() != peptide.length()) {
			// No localization on the representative itself (or a malformed
			// localized sequence that cannot be indexed safely) - fall back to
			// the same no-position bracket the family key already uses.
			return familyKey;
		}
		String bracket = String.format(Locale.ROOT, "[%+.4f]", canonicalMass);
		int insertAt = Math.min(site + 1, peptide.length());
		return peptide.substring(0, insertAt) + bracket + peptide.substring(insertAt);
	}

	/**
	 * One group's final {@code getFlashLfqMonoisotopicMass()} value, from the
	 * same representative PSM and the same delta mass
	 * {@link #decorateWithRepresentativeSite} used to build the "Full
	 * Sequence" bracket - see {@link OpenSearchPsm#getFlashLfqMonoisotopicMass()}'s
	 * javadoc for why this must be {@code calc_neutral_pep_mass} plus the
	 * delta, not {@code calc_neutral_pep_mass} alone, and why it must be
	 * resolved once per group rather than read per PSM: FlashLFQ rejects any
	 * later row whose mass does not exactly match the first one it saw for a
	 * given Full Sequence, and two PSMs of the same peptide+delta almost never
	 * carry bit-identical raw {@code precursor_neutral_mass} values.
	 */
	private static double representativeMonoisotopicMass(OpenSearchPsm representative) {
		if (representative.isClosedStage()) {
			return representative.getCalcNeutralPepMass();
		}
		double canonicalMass = canonicalMassOf(representative);
		if (Double.isNaN(canonicalMass)) {
			// No fitted peak: report the actually-observed ion mass, matching
			// getModificationFamilyKey()'s own [unexplained:+massDiff] branch,
			// which already keys on this PSM's raw massDiff (== precursor mass
			// minus calc_neutral_pep_mass).
			return representative.getPrecursorNeutralMass();
		}
		if (Math.abs(canonicalMass) <= 0.05) {
			return representative.getCalcNeutralPepMass();
		}
		return representative.getCalcNeutralPepMass() + canonicalMass;
	}

	/**
	 * The rematch step: overlays FlashLFQ's real quantification onto the
	 * identification-truth peptides {@link #aggregatePeptides} built.
	 * <p>
	 * {@code validated_psms.tsv} is written once per PSM, per sample, with that
	 * sample's own round-2 database's protein list - not reconciled across
	 * samples (see the comment at {@code MagOpenSearchHeadlessRunner}'s
	 * {@code writeValidatedPsms} call site). FlashLFQ is invoked once across
	 * every sample's rows together, so its own {@code Protein Groups} column
	 * for a peptide reflects every sample's evidence for that peptide, not just
	 * whichever sample this aggregation's representative PSM happened to come
	 * from. Unioning rather than replacing means this can only add proteins
	 * FlashLFQ saw and this aggregation's own psms-derived set did not -
	 * exactly the cross-sample gap that motivated this method - and never lose
	 * one a single PSM's own {@code allProteins} already had.
	 * <p>
	 * Joined on {@link MetaPeptide#getModSeq()}, which is guaranteed to be the
	 * same canonical modified-sequence string on both sides: it is what
	 * {@link #aggregatePeptides} grouped PSMs by, what
	 * {@code MagOpenSearchHeadlessRunner} wrote as the "Full Sequence" column
	 * FlashLFQ read, and what FlashLFQ echoes back as its own "Sequence"
	 * column (see {@link OpenSearchPsm#getCanonicalModifiedSequence()}).
	 */
	private static void rematchWithFlashLfq(MetaPeptide[] peptides,
			Map<String, FlashLfqQuanPeptide> quantifiedByModSeq, String label, ProgressCallback cb) {
		int matched = 0;
		int proteinSetExpanded = 0;
		for (MetaPeptide peptide : peptides) {
			if (peptide == null || peptide.getModSeq() == null) {
				continue;
			}
			FlashLfqQuanPeptide quantified = quantifiedByModSeq.get(peptide.getModSeq());
			if (quantified == null) {
				// Should only happen for a peptide whose only supporting PSM was,
				// for some reason, absent from the idt file FlashLFQ actually read
				// - every validated target PSM aggregatePeptides saw was written
				// to validated_psms.tsv unchanged, so this is a signal something
				// upstream diverged, not an expected steady-state outcome.
				continue;
			}
			matched++;
			peptide.setIntensity(quantified.getIntensity());
			peptide.setIdenType(quantified.getIdenType());

			LinkedHashSet<String> union = new LinkedHashSet<>();
			String[] existing = peptide.getProteins();
			if (existing != null) {
				for (String protein : existing) {
					if (protein != null && !protein.isBlank()) {
						union.add(protein.trim());
					}
				}
			}
			int before = union.size();
			for (String protein : quantified.getProteins()) {
				if (protein != null && !protein.isBlank()) {
					union.add(protein.trim());
				}
			}
			if (union.size() > before) {
				proteinSetExpanded++;
			}
			peptide.setProteins(union.toArray(new String[0]));
		}
		cb.onMessage(label + " - " + matched + "/" + peptides.length
				+ " peptides matched a quantified row (" + (peptides.length - matched)
				+ " kept phase-1 zero intensity), " + proteinSetExpanded
				+ " had their protein set expanded by FlashLFQ's cross-sample evidence");
	}

	/**
	 * Reads FlashLFQ's {@code QuantifiedPeptides.tsv}, keyed by its "Sequence"
	 * column. {@code validated_psms.tsv} keeps decoys for FDR (see
	 * {@link #targetsOnly}), so this holds target and decoy rows alike; FlashLFQ
	 * keys rows by that same full sequence, so a target and a decoy peptide never
	 * share one. {@link FlashLfqQuanPeptide#isDecoy()} cannot tell them apart
	 * here: it tests pFind's {@code REV_}, not {@link TargetDecoyFasta#DECOY_TAG}.
	 */
	static Map<String, FlashLfqQuanPeptide> readQuantifiedPeptides(File quantifiedPeptidesTsv,
			List<String> fileList) {
		FlashLfqQuanPepReader reader = new FlashLfqQuanPepReader(quantifiedPeptidesTsv, new ArrayList<>(fileList));
		Map<String, FlashLfqQuanPeptide> byModSeq = new HashMap<>();
		for (FlashLfqQuanPeptide quantified : reader.getMetaPeptides()) {
			if (quantified == null || quantified.getModSeq() == null) {
				continue;
			}
			byModSeq.put(quantified.getModSeq(), quantified);
		}
		return byModSeq;
	}

	/**
	 * Package-visible (not private): {@link MagOpenSearchHeadlessRunner} reuses this
	 * exact ranking when it writes the "Score" column of the FlashLFQ-format
	 * {@code validated_psms.tsv}, so the score FlashLFQ sees always matches the score
	 * peptide aggregation here picked the representative PSM by.
	 */
	static double score(OpenSearchPsm psm) {
		if (!Double.isNaN(psm.getClassScore())) {
			return psm.getClassScore();
		}
		return Double.isNaN(psm.getHyperscore()) ? 0.0 : psm.getHyperscore();
	}

	/**
	 * Resolves a source file name to its intensity-column index. Matches
	 * {@code DdaMagPostAnalysisService.FileIndexResolver}: MSFragger's pepXML
	 * {@code base_name} is already the sample base name, so a direct lookup is
	 * enough — but fall back to a suffix match in case an extension survives.
	 */
	private static int resolveFileIndex(String fileName, Map<String, Integer> fileIdMap, List<String> fileList) {
		if (fileName == null || fileName.isBlank()) {
			return -1;
		}
		Integer id = fileIdMap.get(fileName);
		if (id != null) {
			return id;
		}
		String stripped = fileName;
		int dot = stripped.lastIndexOf('.');
		if (dot > 0) {
			stripped = stripped.substring(0, dot);
			id = fileIdMap.get(stripped);
			if (id != null) {
				return id;
			}
		}
		for (int i = 0; i < fileList.size(); i++) {
			String candidate = fileList.get(i);
			if (candidate != null && (stripped.endsWith(candidate) || candidate.endsWith(stripped))) {
				return i;
			}
		}
		return -1;
	}

	private static HashMap<String, String> buildRazorProteinMap(MetaPeptide[] peptides, String[] razorProteins) {
		HashMap<String, String> map = new HashMap<>();
		if (peptides == null || razorProteins == null) {
			return map;
		}
		int n = Math.min(peptides.length, razorProteins.length);
		for (int i = 0; i < n; i++) {
			if (peptides[i] != null && razorProteins[i] != null) {
				map.put(peptides[i].getModSeq(), razorProteins[i]);
			}
		}
		return map;
	}
}
