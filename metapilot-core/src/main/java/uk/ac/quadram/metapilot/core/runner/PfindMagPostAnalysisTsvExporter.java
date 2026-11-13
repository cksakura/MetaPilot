package uk.ac.quadram.metapilot.core.runner;

import uk.ac.quadram.metapilot.core.mag.DdaMagPostAnalysisService;
import uk.ac.quadram.metapilot.core.workflow.ProgressCallback;
import uk.ac.quadram.metapilot.dbSearch.pfind.PFindIsobaricTask;
import uk.ac.quadram.metapilot.dbSearch.pfind.PFindResultReader;
import uk.ac.quadram.metapilot.quant.flashLFQ.FlashLfqPsm;
import uk.ac.quadram.metapilot.quant.flashLFQ.FlashLfqPsmReader;
import uk.ac.quadram.metapilot.task.io.qpx.QpxPsmExport;
import uk.ac.quadram.metapilot.quant.flashLFQ.FlashLfqQuanPepReader;
import uk.ac.quadram.metapilot.quant.flashLFQ.FlashLfqQuanPeptide;
import uk.ac.quadram.metapilot.task.io.pep.MetaPeptide;
import uk.ac.quadram.metapilot.task.io.pro.MetaProtein;
import uk.ac.quadram.metapilot.task.mag.MagDbItem;
import uk.ac.quadram.metapilot.task.par.MetaConstants;
import uk.ac.quadram.metapilot.task.par.IsobaricLabelSettings;
import uk.ac.quadram.metapilot.task.par.MetaData;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;

/**
 * Headless post-analysis and TSV export for a DDA MAG pFind workflow, from an
 * already-produced peptide quantification (FlashLFQ label-free, or the
 * isobaric table {@code PFindIsobaricTask} writes). Identification and
 * quantification come from the Windows GUI or
 * {@link PfindMagSearchHeadlessRunner}.
 *
 * Scope/inputs: this expects the same on-disk shape the Windows GUI
 * ({@code DDAMagPfindTask}) produces under {@code <result>/pfind_MAG_search/combined/}:
 * <ul>
 * <li>{@code QuantifiedPeptides.tsv} - FlashLFQ's native label-free peptide output</li>
 * <li>{@code psms.tsv} - the combined pFind PSMs written by {@code OpenPsmTsvWriter}
 * that were fed into FlashLFQ</li>
 * </ul>
 * Only these two files are required to exist; no other intermediate files
 * (marker/sample_db search results, ExperimentalDesign.tsv, etc.) are read or
 * validated - this is intentionally a "resume from quant output" entry point,
 * not a re-run of upstream steps.
 *
 * This mirrors {@link FragpipePostAnalysisTsvExporter}/{@link DiaMagPostAnalysisTsvExporter}
 * in shape and reuses the same shared, platform-independent core services
 * ({@link DdaMagPostAnalysisService}, {@link DdaMagReportExportService}) that back the
 * GUI's {@code MetaPilotMagTask.DDAMagPostAnalysisTask} for the pFind workflow - no
 * JavaFX dependencies.
 */
public final class PfindMagPostAnalysisTsvExporter {

	private PfindMagPostAnalysisTsvExporter() {
	}

	public static final class Result {
		public File idTsv;
		public File finalSummaryTsv;
		public File finalPeptidesTsv;
		public File razorProteinsTsv;
		public File allProteinsTsv;
		public String[] intensityTitles;
		public MetaProtein[] metaProteins;
		public HashMap<String, double[]> genomeIntensityMap;
		public HashMap<String, double[]> genomeUniqueIntensityMap;
		public HashMap<String, Double> genomeScoreMap;
		/** genome -> double[]{S, p_emp, q_value} when PSM-score/top-k scoring is used; null otherwise */
		public HashMap<String, double[]> genomeDiaStats;
		public HashMap<String, Integer> psmCountMap;
		public HashMap<String, HashSet<String>> totalPepCountMap;
		public HashMap<String, HashSet<String>> razorPepCountMap;
		public HashMap<String, HashSet<String>> uniquePepCountMap;
		public HashMap<String, HashSet<String>> totalPepWithILMap;
		public HashMap<String, HashSet<String>> razorPepWithILMap;
		public HashMap<String, HashSet<String>> uniquePepWithILMap;
	}

	public static Result exportFinalTsvs(MagDbItem magDb, File quantifiedPeptidesTsv, File psmsTsv, File outDir,
			MetaData metadata, String taskName, ProgressCallback cb) throws IOException {
		return exportFinalTsvs(magDb, quantifiedPeptidesTsv, psmsTsv, outDir, metadata, taskName, cb,
				MetaConstants.labelFree, false);
	}

	/**
	 * As above, for either quantification mode the pFind workflow produces.
	 * <p>
	 * Label-free reads FlashLFQ's {@code QuantifiedPeptides.tsv} with one
	 * intensity column per raw file. Isobaric reads the peptide table
	 * {@link PFindIsobaricTask#runOnlyPeptide()} writes, with the intensity
	 * columns {@link PFindIsobaricTask#intensityColumnSuffixes} names (raw file
	 * x channel, or sample name when labels are combined) - the same reader and
	 * columns the GUI's {@code DDAMagPfindTask.exportReport()} uses.
	 *
	 * @param quanMode     {@link MetaConstants#labelFree} or
	 *                     {@link MetaConstants#isobaricLabel}; GUI radio texts
	 *                     ("Label-free", "Isobaric") are accepted
	 * @param combineLabel isobaric only: intensities were combined by sample name
	 * @throws IOException when none of the expected intensity columns are found,
	 *                     or every peptide has zero intensity - the two checks the
	 *                     GUI stops on before post-analysis
	 */
	public static Result exportFinalTsvs(MagDbItem magDb, File quantifiedPeptidesTsv, File psmsTsv, File outDir,
			MetaData metadata, String taskName, ProgressCallback cb, String quanMode, boolean combineLabel)
			throws IOException {
		return exportFinalTsvs(magDb, quantifiedPeptidesTsv, psmsTsv, outDir, metadata, taskName, cb, quanMode,
				combineLabel, null);
	}

	/**
	 * As above, with the spectra folder holding {@code ms2Count.txt}, so the
	 * MS/MS column of {@code final_summary.tsv} is the number of MS/MS spectra
	 * and not the number of identified ones, as in the GUI.
	 */
	public static Result exportFinalTsvs(MagDbItem magDb, File quantifiedPeptidesTsv, File psmsTsv, File outDir,
			MetaData metadata, String taskName, ProgressCallback cb, String quanMode, boolean combineLabel,
			File spectraDir) throws IOException {
		if (cb == null)
			cb = new ProgressCallback.NoOp();
		if (magDb == null)
			throw new IllegalArgumentException("magDb is null");
		if (quantifiedPeptidesTsv == null || !quantifiedPeptidesTsv.exists() || quantifiedPeptidesTsv.length() == 0) {
			throw new IOException("QuantifiedPeptides.tsv not found or empty: "
					+ (quantifiedPeptidesTsv == null ? "null" : quantifiedPeptidesTsv.getAbsolutePath()));
		}
		if (psmsTsv == null || !psmsTsv.exists() || psmsTsv.length() == 0) {
			throw new IOException(
					"psms.tsv not found or empty: " + (psmsTsv == null ? "null" : psmsTsv.getAbsolutePath()));
		}
		if (outDir == null)
			throw new IllegalArgumentException("outDir is null");
		if (metadata == null || metadata.getRawFiles() == null || metadata.getRawFiles().length == 0) {
			throw new IllegalArgumentException("metadata with raw files is required to align quantification columns");
		}
		outDir.mkdirs();

		// Build fileList/fileIdMap from metadata raw files. FlashLFQ's own
		// "Intensity_<name>" column headers are keyed by raw file base name (matching
		// the sample folder names the GUI uses), and FlashLfqQuanPepReader already
		// falls back to base-name/path-suffix matching - so this does not depend on
		// ExperimentalDesign.tsv or any other intermediate file existing.
		String[] rawFiles = metadata.getRawFiles();
		ArrayList<String> fileList = new ArrayList<>();
		HashMap<String, Integer> fileIdMap = new HashMap<>();
		for (String rawFile : rawFiles) {
			String base = baseNameNoExt(rawFile);
			fileIdMap.put(base, fileIdMap.size());
			fileList.add(base);
		}
		String mode = IsobaricLabelSettings.normalizeQuantitativeMode(quanMode);
		boolean isobaric = MetaConstants.isobaricLabel.equals(mode);
		if (!isobaric && !MetaConstants.labelFree.equals(mode)) {
			throw new IllegalArgumentException("Unsupported quantitative mode for the pFind workflow: " + quanMode
					+ " (expected " + MetaConstants.labelFree + " or " + MetaConstants.isobaricLabel + ")");
		}
		String[] intensityTitle = isobaric
				? PFindIsobaricTask.intensityColumnSuffixes(metadata, fileList, combineLabel)
				: fileList.toArray(new String[0]);

		cb.onMessage(taskName + ": loading PSMs from " + psmsTsv.getName() + " ...");
		FlashLfqPsmReader psmReader = new FlashLfqPsmReader(psmsTsv);
		FlashLfqPsm[] psms = psmReader.getPSMs();
		if (psms == null)
			psms = new FlashLfqPsm[0];

		cb.onMessage(taskName + ": loading quantified peptides from " + quantifiedPeptidesTsv.getName() + " ...");
		FlashLfqQuanPepReader quanPepReader = isobaric
				? FlashLfqQuanPepReader.withIntensityColumns(quantifiedPeptidesTsv, fileList, intensityTitle)
				: new FlashLfqQuanPepReader(quantifiedPeptidesTsv, fileList);
		cb.onMessage(taskName + ": peptide input (" + mode + "): " + quanPepReader.getIntensityMatchSummary());
		if (quanPepReader.getMatchedIntensityCount() == 0) {
			throw new IOException("None of the " + quanPepReader.getExpectedIntensityCount()
					+ " expected intensity columns were found in " + quantifiedPeptidesTsv.getAbsolutePath()
					+ "; peptide quantification cannot be transferred. " + quanPepReader.getIntensityMatchSummary());
		}
		if (quanPepReader.getMatchedIntensityCount() < quanPepReader.getExpectedIntensityCount()) {
			ArrayList<String> missing = quanPepReader.getUnmatchedIntensityTitles();
			cb.onMessage(taskName + ": WARNING: only " + quanPepReader.getMatchedIntensityCount() + "/"
					+ quanPepReader.getExpectedIntensityCount()
					+ " expected intensity columns were found; unmatched columns will be zero. First missing: "
					+ missing.subList(0, Math.min(5, missing.size())));
		}
		FlashLfqQuanPeptide[] quanPeps = quanPepReader.getMetaPeptides();
		if (quanPeps == null)
			quanPeps = new FlashLfqQuanPeptide[0];
		if (quanPeps.length > 0 && quanPepReader.getPeptidesWithIntensity() == 0) {
			throw new IOException("All " + quanPeps.length + " peptides read from "
					+ quantifiedPeptidesTsv.getAbsolutePath() + " have zero intensity; stopping before post-analysis. "
					+ quanPepReader.getIntensityMatchSummary());
		}

		cb.onMessage(taskName + ": peptides=" + quanPeps.length + ", psms=" + psms.length);

		DdaMagPostAnalysisService.Result analysis = DdaMagPostAnalysisService.analyze(
				new DdaMagPostAnalysisService.Input(quanPeps, psms, fileIdMap, fileList, PFindResultReader.REV,
						magDb.getIdentifier()));
        analysis.finalSelection.write(outDir);
        uk.ac.quadram.metapilot.core.mag.MagResultFinalisationService.writeCatalogueContext(outDir, magDb);

		HashMap<String, String> razorProteinMap = buildRazorProteinMap(quanPeps, analysis.razorProteins);
		cb.onMessage(taskName + ": genome scoring="
				+ (analysis.genomeScoring.diaStyle ? "PSM-score/top-k/q-value" : "legacy fallback") + ", genomes="
				+ analysis.genomeScoreMap.size());

		File combinedDir = new File(outDir, "combined");
		File idOut = new File(combinedDir, "id.tsv");
		File summaryOut = new File(outDir, "final_summary.tsv");
		File pepOut = new File(outDir, "final_peptides.tsv");
		File proOut = new File(outDir, "razor_proteins.tsv");
		File allProOut = new File(outDir, "all_proteins.tsv");

		cb.onMessage(taskName + ": writing " + idOut.getName() + " / " + summaryOut.getName());
		DdaMagReportExportService.writeIdentificationSummary(combinedDir, outDir, spectraDir, fileIdMap,
				analysis.filePsmCount, analysis.fileSequenceSets, null);

		cb.onMessage(taskName + ": writing " + pepOut.getName());
		DdaMagReportExportService.writePeptideReport(pepOut, analysis.peptideList, analysis.psmMap,
				analysis.pepPsmCountMap, analysis.pepChargeMap, analysis.proteinIdMap, fileList,
				mode, intensityTitle, razorProteinMap);

		cb.onMessage(taskName + ": writing " + proOut.getName());
		DdaMagReportExportService.writeProteinReport(proOut, analysis.proteins, analysis.proteinIntensityMap,
				analysis.proteinTotalIntenMap, analysis.totalPepCountMap, analysis.razorPepCountMap,
				analysis.psmCountMap, mode, intensityTitle, analysis.peptideList,
				analysis.allProteinScoreMap, analysis.psmMap);

		cb.onMessage(taskName + ": writing " + allProOut.getName());
		DdaMagReportExportService.writeAllProteinReport(allProOut, analysis.allProteins,
				analysis.allProteinIntensityMap, analysis.allProteinTotalIntenMap, analysis.allProteinPepCountMap,
				analysis.razorPepCountMap, analysis.allProteinPsmCountMap, mode, intensityTitle,
				analysis.peptideList, analysis.allProteinScoreMap, analysis.psmMap);

		// The PSM layer as QPX psm.parquet, beside final_peptides.tsv (never fails the export).
		QpxPsmExport.writeQuietly(outDir, psms, cb);

		Result res = new Result();
		res.idTsv = idOut;
		res.finalSummaryTsv = summaryOut;
		res.finalPeptidesTsv = pepOut;
		res.razorProteinsTsv = proOut;
		res.allProteinsTsv = allProOut;
		res.intensityTitles = intensityTitle;
		res.metaProteins = analysis.metaProteins;
		res.genomeIntensityMap = analysis.genomeIntensityMap;
		res.genomeUniqueIntensityMap = analysis.genomeUniqueIntensityMap;
		res.genomeScoreMap = analysis.genomeScoreMap;
		res.genomeDiaStats = analysis.genomeDiaStats;
		res.psmCountMap = analysis.psmCountMap;
		res.totalPepCountMap = analysis.totalPepCountMap;
		res.razorPepCountMap = analysis.razorPepCountMap;
		res.uniquePepCountMap = analysis.uniquePepCountMap;
		res.totalPepWithILMap = analysis.totalPepWithILMap;
		res.razorPepWithILMap = analysis.razorPepWithILMap;
		res.uniquePepWithILMap = analysis.uniquePepWithILMap;
		return res;
	}

	static HashMap<String, String> buildRazorProteinMap(MetaPeptide[] peptides, String[] razorProteins) {
		HashMap<String, String> map = new HashMap<>();
		if (peptides == null || razorProteins == null) {
			return map;
		}
		for (int i = 0; i < peptides.length && i < razorProteins.length; i++) {
			if (peptides[i] == null || razorProteins[i] == null) {
				continue;
			}
			if (peptides[i].getModSeq() != null) {
				map.put(peptides[i].getModSeq(), razorProteins[i]);
			}
			// Exact peptidoform keys only: a modified peptide's base sequence may
			// equal another peptide's full key and overwrite its razor assignment.
		}
		return map;
	}

	private static String baseNameNoExt(String path) {
		if (path == null)
			return "sample";
		String name = new File(path).getName();
		int dot = name.lastIndexOf('.');
		return dot > 0 ? name.substring(0, dot) : name;
	}
}
