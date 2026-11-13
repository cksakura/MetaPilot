package uk.ac.quadram.metapilot.core.runner;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.lang.management.ManagementFactory;
import java.nio.file.Files;
import java.nio.file.StandardCopyOption;
import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.text.ParseException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Locale;

import uk.ac.quadram.metapilot.core.mag.MagGenomeSelectionService;
import uk.ac.quadram.metapilot.core.mag.MagResultLayoutUtil;
import uk.ac.quadram.metapilot.core.resources.MetaPilotResourceConfig;
import uk.ac.quadram.metapilot.core.variant.PFindVariantCandidateSource;
import uk.ac.quadram.metapilot.core.variant.PangenomeVariantCheck;
import uk.ac.quadram.metapilot.core.variant.PangenomeVariantSettings;
import uk.ac.quadram.metapilot.core.variant.VariantCandidate;
import uk.ac.quadram.metapilot.core.workflow.ProgressCallback;
import uk.ac.quadram.metapilot.dbSearch.open.io.OpenPsmTsvWriter;
import uk.ac.quadram.metapilot.dbSearch.pfind.PFindIsobaricTask;
import uk.ac.quadram.metapilot.dbSearch.pfind.PFindPSM;
import uk.ac.quadram.metapilot.dbSearch.pfind.PFindPsmReader;
import uk.ac.quadram.metapilot.dbSearch.pfind.PFindResultReader;
import uk.ac.quadram.metapilot.dbSearch.pfind.PFindTask;
import uk.ac.quadram.metapilot.core.mod.IsobaricTag;
import uk.ac.quadram.metapilot.quant.flashLFQ.FlashLfqTask;
import uk.ac.quadram.metapilot.task.dda.mag.par.DDAMagPfindParameters;
import uk.ac.quadram.metapilot.task.io.par.MetaOpenPar;
import uk.ac.quadram.metapilot.task.mag.MagDbItem;
import uk.ac.quadram.metapilot.task.par.IsobaricLabelSettings;
import uk.ac.quadram.metapilot.task.par.MetaConstants;
import uk.ac.quadram.metapilot.task.par.MetaData;

/**
 * Headless, end-to-end DDA MAG pFind workflow - the command-line counterpart of
 * the Windows GUI's {@code DDAMagPfindTask}, with no JavaFX dependency.
 * <p>
 * <b>Windows only.</b> pFind ({@code Searcher.exe}, {@code pParse}/
 * {@code pParse2Plus}) and FlashLFQ's {@code CMD.exe} are Windows executables,
 * and {@link PFindTask} drives them through {@code .bat} files. {@link #run}
 * refuses to start on any other operating system.
 * <p>
 * Steps, in the GUI's order, each skipped when its output already exists so an
 * interrupted run resumes where it stopped:
 * <ol>
 * <li>Spectra conversion: raw to {@code .pf2}/{@code .ms2} with pFind's
 * converter, retrying failed files without pParse2Plus; then {@code rt.txt}
 * and {@code ms2Count.txt} from the {@code .ms2} files.</li>
 * <li>Marker search: each raw file against the MAG catalogue's marker
 * FASTA.</li>
 * <li>Genome selection per raw file and its sample-specific FASTA (selected
 * genomes, plus the host database when appended).</li>
 * <li>Sample-specific search: each raw file against its own FASTA.</li>
 * <li>PSM combination into {@code combined/psms.tsv}.</li>
 * <li>Quantification: FlashLFQ label-free; for isobaric labelling, FlashLFQ
 * followed by {@link PFindIsobaricTask#runOnlyPeptide()}.</li>
 * <li>Post-analysis, taxonomy, functional annotation and the HTML report,
 * through {@link PfindMagHeadlessRunner#postAnalysis} - the same code
 * {@code pfind-mag-report} runs.</li>
 * <li>Pangenome variant check, open search only.</li>
 * </ol>
 * The result layout is the GUI's: {@code <result>/pfind_MAG_search/} with
 * per-raw-file folders under {@code search_results/}, each holding
 * {@code marker/} and {@code sample_db/}, plus {@code combined/}; spectra go
 * to {@code <result>/spectra/}. A folder started in the GUI can be finished
 * here and the reverse.
 * <p>
 * Differences from the GUI, all deliberate:
 * <ul>
 * <li>Decoy PSMs are written to {@code psms.tsv} once. The GUI's
 * {@code combineResults()} writes every decoy that has a retention time twice,
 * which doubles the decoys downstream FDR and genome scoring see.</li>
 * <li>MGF files are not written. In the GUI they feed only the DBReducer path,
 * which {@code DDAMagPfindTask} never enables.</li>
 * <li>Only Thermo {@code .raw} input is accepted, and all raw files must share
 * one folder (FlashLFQ reads them from a single folder). The GUI accepts other
 * input and then fails later.</li>
 * <li>Each sample FASTA is written to a temporary file and renamed when
 * complete, so an interrupted write is not reused as a finished database.</li>
 * <li>Post-analysis always runs; the GUI skips it when its outputs are newer
 * than {@code psms.tsv} and the peptide table.</li>
 * </ul>
 */
public final class PfindMagSearchHeadlessRunner {

	public static final String RESULT_DIR = "pfind_MAG_search";
	static final String MARKER_SEARCH_DIR = "marker";
	static final String SAMPLE_DB_SEARCH_DIR = "sample_db";
	static final String FILTERED_SPECTRA = "pFind-Filtered.spectra";
	private static final String TASK = "DDA MAG pFind";

	private PfindMagSearchHeadlessRunner() {
	}

	/** Everything a run needs besides the parameter file. */
	public static final class Request {
		public DDAMagPfindParameters params;
		/** {@code pFind.exe} (resource key PFIND); its folder must hold {@code Searcher.exe}. */
		public File pFindExe;
		/** FlashLFQ {@code CMD.exe} (resource key FLASHLFQ). */
		public File flashLfqExe;
		/** MAG DB registry JSON (MAGDB_1_0.json). */
		public File magdbConfigJson;
		/** {@code func_def.db}; null skips functional annotation. */
		public File funcDbFile;
		/**
		 * Raw files searched at the same time. 0 = one per 5 GB of physical memory,
		 * as the GUI decides it.
		 */
		public int parallelTasks;
		public boolean generateHtmlReport = true;
		public String reportPythonHint;
		public MetaPilotResourceConfig reportResourceConfig;
		public boolean exportMetaAtlasParquet;
		public ProgressCallback cb;
	}

	public static final class Result {
		public final File resultRoot;
		public final File psmsTsv;
		public final File quantifiedPeptidesTsv;
		public final PfindMagPostAnalysisTsvExporter.Result post;

		Result(File resultRoot, File psmsTsv, File quantifiedPeptidesTsv, PfindMagPostAnalysisTsvExporter.Result post) {
			this.resultRoot = resultRoot;
			this.psmsTsv = psmsTsv;
			this.quantifiedPeptidesTsv = quantifiedPeptidesTsv;
			this.post = post;
		}
	}

	/** true on Windows, the only platform pFind runs on. */
	public static boolean isSupportedPlatform() {
		return System.getProperty("os.name", "").toLowerCase(Locale.ROOT).startsWith("windows");
	}

	public static Result run(Request req) throws Exception {
		if (req == null || req.params == null) {
			throw new IllegalArgumentException("params is null");
		}
		ProgressCallback cb = req.cb != null ? req.cb : new ProgressCallback.NoOp();
		if (!isSupportedPlatform()) {
			throw new IllegalStateException("The pFind workflow runs on Windows only (pFind and FlashLFQ are Windows "
					+ "executables); this is " + System.getProperty("os.name") + ". Use pfind-mag-report to "
					+ "post-process results produced on Windows.");
		}
		Run run = new Run(req, cb);
		return run.execute();
	}

	/** One workflow run; holds the layout derived from the parameters. */
	private static final class Run {
		final Request req;
		final DDAMagPfindParameters params;
		final ProgressCallback cb;

		MagDbItem mag;
		File pFindBin;
		File resultRoot;
		File combinedDir;
		File spectraDir;
		File rawDir;
		String[] rawFiles;
		String[] rawNames;
		File[] sampleDirs;
		String[] pf2Files;
		int maxTaskCount;
		String quanMode;

		Run(Request req, ProgressCallback cb) {
			this.req = req;
			this.params = req.params;
			this.cb = cb;
		}

		Result execute() throws Exception {
			prepare();

			step("Step 1/8: spectra conversion");
			convertSpectra();

			step("Step 2/8: marker database search");
			searchMarker();

			step("Step 3/8: genome selection and sample-specific databases");
			buildSampleDatabases();

			step("Step 4/8: sample-specific database search");
			searchSampleDatabases();

			step("Step 5/8: combining PSMs");
			File psms = combinePsms();

			step("Step 6/8: quantification (" + quanMode + ")");
			File quanPep = quantify(psms);

			step("Step 7/8: post-analysis, taxonomy, functional annotation and report");
			PfindMagPostAnalysisTsvExporter.Result post = PfindMagHeadlessRunner.postAnalysis(params, mag, resultRoot,
					quanPep, psms, req.funcDbFile, cb, req.generateHtmlReport, req.reportPythonHint,
					req.reportResourceConfig, req.exportMetaAtlasParquet);

			if (params.isOpenSearch()) {
				step("Step 8/8: pangenome variant check");
				runPangenomeVariantCheck();
			} else {
				info("Step 8/8: pangenome variant check runs only for open search; not run for this closed search");
			}

			info("completed: " + resultRoot.getAbsolutePath());
			cb.onComplete(resultRoot.getAbsolutePath());
			return new Result(resultRoot, psms, quanPep, post);
		}

		// ---------------------------------------------------------------- setup

		void prepare() throws IOException {
			if (req.pFindExe == null || !req.pFindExe.isFile()) {
				throw new IOException("pFind executable not found: " + describe(req.pFindExe));
			}
			pFindBin = req.pFindExe.getParentFile();
			File searcher = new File(pFindBin, "Searcher.exe");
			File templateParam = new File(pFindBin, "template" + File.separator + "param");
			for (File required : new File[] { searcher, new File(templateParam, "pFind.cfg"),
					new File(templateParam, "pParse.cfg"), new File(templateParam, "pQuant.cfg") }) {
				if (!required.isFile()) {
					throw new IOException("pFind installation is incomplete, missing: " + required.getAbsolutePath());
				}
			}
			if (req.flashLfqExe == null || !req.flashLfqExe.isFile()) {
				throw new IOException("FlashLFQ executable not found: " + describe(req.flashLfqExe));
			}

			quanMode = IsobaricLabelSettings
					.normalizeQuantitativeMode(params.getQuantitativeMode());
			if (quanMode == null || quanMode.isBlank()) {
				quanMode = MetaConstants.labelFree;
			}
			if (!MetaConstants.labelFree.equals(quanMode) && !MetaConstants.isobaricLabel.equals(quanMode)) {
				throw new IllegalArgumentException("The pFind workflow supports label-free and isobaric "
						+ "quantification; the parameter file asks for: " + params.getQuantitativeMode());
			}

			mag = PfindMagHeadlessRunner.resolveMagDb(params.getMagDbCatalogueId(), params.getMagDbVersion(),
					req.magdbConfigJson);
			if (mag == null) {
				throw new IllegalArgumentException("MAG database not found in " + describe(req.magdbConfigJson)
						+ " for catalogue=" + params.getMagDbCatalogueId() + " version=" + params.getMagDbVersion());
			}
			File markerFasta = mag.getHapFastaFile(params.getGenomeSelectorMode());
			if (markerFasta == null || !markerFasta.isFile()) {
				throw new IOException("MAG database marker FASTA not found for genome selector "
						+ params.getGenomeSelectorMode() + ": " + describe(markerFasta));
			}
			if (mag.getCurrentFile() == null || !mag.getCurrentFile().exists()) {
				throw new IOException("MAG database folder not found: " + describe(mag.getCurrentFile()));
			}
			params.setMagDbItem(mag);
			info("MAG DB: " + mag.getCatalogueID() + " / " + mag.getUsedVersion().getCatalogueVersion());

			if (params.isAppendHostDb()) {
				String host = params.getHostDB();
				if (host == null || host.isBlank() || !new File(host).isFile()) {
					throw new IOException("append_host_db is set but the host database was not found: " + host);
				}
			}

			MetaData metadata = params.getMetadata();
			rawFiles = metadata == null ? null : metadata.getRawFiles();
			if (rawFiles == null || rawFiles.length == 0) {
				throw new IllegalArgumentException("The parameter file lists no raw files");
			}
			if (params.getResult() == null || params.getResult().isBlank()) {
				throw new IllegalArgumentException("The parameter file has no output_directory");
			}
			rawNames = new String[rawFiles.length];
			HashSet<String> seen = new HashSet<>();
			for (int i = 0; i < rawFiles.length; i++) {
				File raw = new File(rawFiles[i]);
				String name = raw.getName();
				if (!name.toLowerCase(Locale.ROOT).endsWith(".raw")) {
					throw new IllegalArgumentException("The pFind workflow takes Thermo .raw files: " + rawFiles[i]);
				}
				if (!raw.isFile()) {
					throw new IOException("Raw file not found: " + raw.getAbsolutePath());
				}
				File parent = raw.getAbsoluteFile().getParentFile();
				if (rawDir == null) {
					rawDir = parent;
				} else if (!rawDir.equals(parent)) {
					throw new IllegalArgumentException("All raw files must be in one folder (FlashLFQ reads them "
							+ "from a single folder): " + rawDir + " and " + parent);
				}
				rawNames[i] = name.substring(0, name.length() - ".raw".length());
				if (!seen.add(rawNames[i].toLowerCase(Locale.ROOT))) {
					throw new IllegalArgumentException("Raw file name listed twice: " + rawNames[i]);
				}
			}

			if (MetaConstants.isobaricLabel.equals(quanMode)) {
				prepareIsobaricLabels(metadata);
			}

			resultRoot = new File(params.getResult(), RESULT_DIR);
			resultRoot.mkdirs();
			params.setDbSearchResultFile(resultRoot);
			combinedDir = new File(resultRoot, "combined");
			combinedDir.mkdirs();
			File samplesRoot = MagResultLayoutUtil.resolveSearchResultsRoot(resultRoot);
			sampleDirs = new File[rawFiles.length];
			for (int i = 0; i < rawFiles.length; i++) {
				sampleDirs[i] = new File(samplesRoot, rawNames[i]);
				new File(sampleDirs[i], MARKER_SEARCH_DIR).mkdirs();
			}
			spectraDir = params.getSpectraFile();

			String suffix = ms2Suffix(params.getMs2ScanMode());
			pf2Files = new String[rawFiles.length];
			for (int i = 0; i < rawFiles.length; i++) {
				pf2Files[i] = new File(spectraDir, rawNames[i] + suffix + ".pf2").getAbsolutePath();
			}

			maxTaskCount = req.parallelTasks > 0 ? req.parallelTasks : memoryBasedTaskCount();

			info("result folder: " + resultRoot.getAbsolutePath());
			info("raw files: " + rawFiles.length + " in " + rawDir.getAbsolutePath());
			info("MS2 scan mode: " + params.getMs2ScanMode() + ", search: "
					+ (params.isOpenSearch() ? "open" : "closed") + ", enzyme: " + params.getEnzyme()
					+ ", missed cleavages: " + params.getMissCleavages() + ", digestion: "
					+ params.getDigestModeString());
			info("quantification: " + quanMode + ", threads: " + params.getThreadCount()
					+ ", parallel searches: " + maxTaskCount);
			info("pFind: " + req.pFindExe.getAbsolutePath() + ", FlashLFQ: " + req.flashLfqExe.getAbsolutePath());
		}

		/**
		 * The channel labels isobaric quantification indexes by. The GUI sets them
		 * from its isobaric table when it builds the run; a parameter file restores
		 * the table but not the labels, so they are rebuilt from it here, falling
		 * back to the resolved tag's channels.
		 */
		void prepareIsobaricLabels(MetaData metadata) {
			if (metadata.getLabelCount() > 0) {
				return;
			}
			ArrayList<String> labels = new ArrayList<>();
			List<String[]> table = params.getIsobaricTableData();
			if (table != null) {
				for (String[] row : table) {
					if (row != null && row.length > 0 && row[0] != null && !row[0].isBlank()) {
						labels.add(row[0].trim());
					}
				}
			}
			if (labels.isEmpty()) {
				IsobaricTag tag = params.getIsobaricTag();
				if (tag != null) {
					labels.addAll(Arrays.asList(tag.getReportTitle()));
				}
			}
			if (labels.isEmpty()) {
				throw new IllegalArgumentException("Isobaric quantification is selected but the parameter file has "
						+ "no isobaric type or channel table");
			}
			metadata.setLabelTitle(labels.toArray(new String[0]));
			info("isobaric channels: " + labels.size() + " (" + String.join(", ", labels) + ")");
		}

		// ------------------------------------------------------- 1. spectra

		void convertSpectra() throws IOException {
			ArrayList<String> toExtract = new ArrayList<>();
			for (int i = 0; i < rawFiles.length; i++) {
				if (!hasContent(new File(pf2Files[i]))) {
					toExtract.add(rawFiles[i]);
				}
			}
			if (toExtract.isEmpty()) {
				info("spectra already converted in " + spectraDir.getAbsolutePath());
			} else {
				info("converting " + toExtract.size() + " raw file(s) to " + spectraDir.getAbsolutePath());
				PFindTask converter = new PFindTask(pFindBin, Math.max(1, params.getThreadCount()));
				converter.extract(spectraDir.getAbsolutePath(), toExtract.toArray(new String[0]));
				converter.run(toExtract.size());

				ArrayList<String> retry = new ArrayList<>();
				for (String raw : toExtract) {
					if (!hasContent(new File(pf2Files[indexOf(raw)]))) {
						retry.add(raw);
					}
				}
				if (!retry.isEmpty()) {
					info("no .pf2 for " + retry.size() + " file(s); retrying without pParse2Plus");
					converter.extractWithoutPParse2Plus(spectraDir.getAbsolutePath(), retry.toArray(new String[0]));
					converter.run(retry.size());
				}

				ArrayList<String> missing = new ArrayList<>();
				for (int i = 0; i < rawFiles.length; i++) {
					if (!hasContent(new File(pf2Files[i]))) {
						missing.add(rawFiles[i] + " -> expected " + pf2Files[i]);
					}
				}
				if (!missing.isEmpty()) {
					throw new IOException(missing.size() + " raw file(s) have no .pf2 spectra after conversion and "
							+ "retry (check MS2 scan mode " + params.getMs2ScanMode() + "):\n  "
							+ String.join("\n  ", missing));
				}
			}
			writeRetentionTimes();
		}

		/** rt.txt (minutes) and ms2Count.txt from the .ms2 files, as the GUI writes them. */
		void writeRetentionTimes() throws IOException {
			File rtFile = new File(spectraDir, "rt.txt");
			File ms2CountFile = new File(spectraDir, "ms2Count.txt");
			if (hasContent(rtFile) && hasContent(ms2CountFile)) {
				return;
			}
			info("writing retention times to " + rtFile.getAbsolutePath());
			File rtTmp = new File(spectraDir, "rt.txt.tmp");
			File countTmp = new File(spectraDir, "ms2Count.txt.tmp");
			try (PrintWriter rtWriter = new PrintWriter(rtTmp); PrintWriter countWriter = new PrintWriter(countTmp)) {
				rtWriter.println("FileName\tScannum\tRt");
				countWriter.println("FileName\tMS2 count");
				DecimalFormat df = null;
				for (String rawName : rawNames) {
					File ms2 = new File(spectraDir, rawName + ".ms2");
					if (!ms2.isFile()) {
						throw new IOException("MS2 file not found after conversion: " + ms2.getAbsolutePath());
					}
					int scan = 0;
					int count = 0;
					try (BufferedReader reader = new BufferedReader(new FileReader(ms2))) {
						String line;
						while ((line = reader.readLine()) != null) {
							String[] cs = line.split("\t");
							if (cs.length < 3) {
								continue;
							}
							if (cs[0].equals("S")) {
								if (df == null && cs.length > 3) {
									df = new DecimalFormat("#,##0.00", new DecimalFormatSymbols(
											cs[3].contains(",") ? Locale.GERMANY : Locale.US));
								}
								scan = Integer.parseInt(cs[1]);
								count++;
							} else if (cs[0].equals("I") && cs[1].equals("RetTime") && df != null) {
								double rtSeconds;
								try {
									rtSeconds = df.parse(cs[2]).doubleValue();
								} catch (ParseException e) {
									throw new IOException("Cannot parse retention time in " + ms2 + ": " + line, e);
								}
								rtWriter.println(rawName + "\t" + scan + "\t" + rtSeconds / 60.0);
							}
						}
					}
					countWriter.println(rawName + "\t" + count);
				}
			}
			Files.move(rtTmp.toPath(), rtFile.toPath(), StandardCopyOption.REPLACE_EXISTING);
			Files.move(countTmp.toPath(), ms2CountFile.toPath(), StandardCopyOption.REPLACE_EXISTING);
		}

		// ------------------------------------------------------- 2. marker search

		void searchMarker() throws IOException {
			File markerFasta = mag.getHapFastaFile(params.getGenomeSelectorMode());
			ArrayList<Integer> todo = new ArrayList<>();
			for (int i = 0; i < rawFiles.length; i++) {
				if (!hasContent(new File(new File(sampleDirs[i], MARKER_SEARCH_DIR), FILTERED_SPECTRA))) {
					todo.add(i);
				}
			}
			if (todo.isEmpty()) {
				info("marker search results already exist for all " + rawFiles.length + " file(s)");
			} else {
				int pool = Math.min(maxTaskCount, todo.size());
				PFindTask task = newSearchTask(pool);
				for (int i : todo) {
					File markerDir = new File(sampleDirs[i], MARKER_SEARCH_DIR);
					markerDir.mkdirs();
					addSearch(task, i, markerDir, markerFasta);
				}
				info("marker search: " + todo.size() + " file(s) against " + markerFasta.getName() + ", " + pool
						+ " at a time");
				task.run(pool, sampleDirs.length * 12);
			}
			for (File sampleDir : sampleDirs) {
				File markerDir = new File(sampleDir, MARKER_SEARCH_DIR);
				PFindTask.cleanPFindResultFolder(markerDir, true);
				if (!new PFindResultReader(markerDir, spectraDir).hasResult()) {
					throw new IOException("pFind marker search produced no result in " + markerDir.getAbsolutePath());
				}
			}
		}

		// ------------------------------------------- 3. sample-specific FASTA

		void buildSampleDatabases() throws Exception {
			boolean marker4 = "Marker4".equalsIgnoreCase(params.getGenomeSelectorMode());
			String prefix = marker4 ? "marker_genome" : "genome";
			for (File sampleDir : sampleDirs) {
				File markerDir = new File(sampleDir, MARKER_SEARCH_DIR);
				File dbFile = sampleDatabase(sampleDir);
				if (hasContent(dbFile)) {
					info("sample-specific database for " + sampleDir.getName() + " already exists");
					continue;
				}
				File spectraFile = new File(markerDir, FILTERED_SPECTRA);
				if (!hasContent(spectraFile)) {
					throw new IOException("pFind marker search result not found: " + spectraFile.getAbsolutePath());
				}
				PFindPSM[] psms = new PFindPsmReader(spectraFile).getPsms();

				MagGenomeSelectionService.Request request = new MagGenomeSelectionService.Request();
				request.magDb = mag;
				request.genomeSelectorMode = params.getGenomeSelectorMode();
				request.genomeSelectionStrategy = params.getGenomeSelectionStrategy();
				request.genomeSelectionPsmCoverage = params.getGenomeSelectionPsmCoverage();
				request.genomeSelectionStepGuard = params.getGenomeSelectionStepGuard();
				request.genomeSelectionMaxGenomes = params.getGenomeSelectionMaxGenomes();
				request.psms = psms == null ? Collections.emptyList() : Arrays.asList(psms);
				request.evidenceFile = new File(markerDir, prefix + "_evidence.tsv");
				request.rankedFile = new File(markerDir, prefix + "_ranked.tsv");
				request.selectedFile = new File(markerDir, prefix + "_selected.tsv");
				request.maximumSelected = MagGenomeSelectionService.DEFAULT_MAXIMUM_SELECTED;
				request.sourceDescription = spectraFile.getAbsolutePath();
				HashSet<String> genomes = MagGenomeSelectionService.selectGenomes(request, cb);
				info(sampleDir.getName() + ": " + genomes.size() + " genomes selected ("
						+ (marker4 ? "Marker4" : "Marker2") + ")");
				if (genomes.isEmpty()) {
					info("WARNING: no genomes were selected for " + sampleDir.getName()
							+ "; it will have no sample-specific search");
				}

				File tmp = new File(dbFile.getParentFile(), dbFile.getName() + ".tmp");
				File dbFolder = mag.getOriginalDbFolder();
				int missing = 0;
				try (PrintWriter writer = new PrintWriter(tmp)) {
					for (String genome : genomes) {
						File genomeDb = new File(dbFolder, genome + ".faa");
						if (!genomeDb.exists()) {
							genomeDb = new File(dbFolder, genome + ".fasta");
						}
						if (!genomeDb.exists()) {
							missing++;
							continue;
						}
						copyLines(genomeDb, writer);
					}
					if (params.isAppendHostDb()) {
						copyLines(new File(params.getHostDB()), writer);
					}
				}
				if (missing > 0) {
					info("WARNING: " + missing + " selected genome(s) have no protein FASTA in "
							+ dbFolder.getAbsolutePath());
				}
				Files.move(tmp.toPath(), dbFile.toPath(), StandardCopyOption.REPLACE_EXISTING);
			}
			try {
				uk.ac.quadram.metapilot.task.mag.GenomeSelectionSummary.combine(sampleDirs, MARKER_SEARCH_DIR,
						new File(resultRoot, uk.ac.quadram.metapilot.task.mag.GenomeSelectionSummary.FILE_NAME));
			} catch (IOException e) {
				info("WARNING: could not write genome_selection_summary.tsv: " + e.getMessage());
			}
		}

		// ------------------------------------------- 4. sample-specific search

		void searchSampleDatabases() throws IOException {
			ArrayList<Integer> todo = new ArrayList<>();
			for (int i = 0; i < sampleDirs.length; i++) {
				File sampleSearchDir = new File(sampleDirs[i], SAMPLE_DB_SEARCH_DIR);
				if (new PFindResultReader(sampleSearchDir, spectraDir).hasResult()) {
					info("sample-specific search result already exists in " + sampleSearchDir.getAbsolutePath());
					continue;
				}
				File db = sampleDatabase(sampleDirs[i]);
				if (!db.exists()) {
					throw new IOException("sample-specific database not found: " + db.getAbsolutePath());
				}
				if (db.length() == 0) {
					info("WARNING: " + db.getName() + " is empty; no peptides can be identified from "
							+ sampleDirs[i].getName());
					continue;
				}
				todo.add(i);
			}
			if (!todo.isEmpty()) {
				int pool = Math.min(maxTaskCount, todo.size());
				PFindTask task = newSearchTask(pool);
				for (int i : todo) {
					File sampleSearchDir = new File(sampleDirs[i], SAMPLE_DB_SEARCH_DIR);
					sampleSearchDir.mkdirs();
					addSearch(task, i, sampleSearchDir, sampleDatabase(sampleDirs[i]));
				}
				info("sample-specific search: " + todo.size() + " file(s), " + pool + " at a time");
				task.run(pool, sampleDirs.length * 12);
			}
			int withResult = 0;
			for (File sampleDir : sampleDirs) {
				File sampleSearchDir = new File(sampleDir, SAMPLE_DB_SEARCH_DIR);
				if (sampleSearchDir.exists()) {
					PFindTask.cleanPFindResultFolder(sampleSearchDir, true);
				}
				PFindTask.cleanPFindResultFolder(new File(sampleDir, MARKER_SEARCH_DIR), true);
				if (hasContent(new File(sampleSearchDir, FILTERED_SPECTRA))) {
					withResult++;
				} else {
					info("WARNING: no sample-specific search result for " + sampleDir.getName());
				}
			}
			if (withResult == 0) {
				throw new IOException("The sample-specific search produced no result for any raw file");
			}
		}

		// ------------------------------------------------------- 5. combine

		File combinePsms() throws IOException {
			File combined = new File(combinedDir, "psms.tsv");
			if (combined.exists() && combined.length() > 200L) {
				info("PSMs already combined in " + combined.getAbsolutePath());
				return combined;
			}
			HashMap<String, Double> rtMap = new HashMap<>();
			File rtFile = new File(spectraDir, "rt.txt");
			try (BufferedReader reader = new BufferedReader(new FileReader(rtFile))) {
				String line = reader.readLine();
				while ((line = reader.readLine()) != null) {
					String[] cs = line.split("\t");
					rtMap.put(cs[0] + "_" + cs[1], Double.parseDouble(cs[2]));
				}
			}

			String microbialPrefix = mag.getIdentifier();
			ArrayList<PFindPSM> targets = new ArrayList<>();
			ArrayList<PFindPSM> hosts = new ArrayList<>();
			ArrayList<PFindPSM> decoys = new ArrayList<>();
			for (File sampleDir : sampleDirs) {
				File spectraFile = new File(new File(sampleDir, SAMPLE_DB_SEARCH_DIR), FILTERED_SPECTRA);
				if (!hasContent(spectraFile)) {
					continue;
				}
				PFindPSM[] psms = new PFindPsmReader(spectraFile).getPsms();
				for (PFindPSM psm : psms) {
					if (psm.getBaseSequence().contains("X")) {
						continue;
					}
					ArrayList<String> rev = new ArrayList<>();
					ArrayList<String> host = new ArrayList<>();
					for (String protein : psm.getProteins()) {
						if (protein.startsWith(PFindResultReader.REV)) {
							rev.add(protein);
						} else if (!protein.startsWith(microbialPrefix)) {
							host.add(protein);
						}
					}
					if (!rev.isEmpty()) {
						psm.setProteins(rev.toArray(new String[0]));
						decoys.add(psm);
					} else if (!host.isEmpty()) {
						psm.setProteins(host.toArray(new String[0]));
						hosts.add(psm);
					} else {
						targets.add(psm);
					}
				}
				info(sampleDir.getName() + ": " + psms.length + " PSMs");
			}
			if (targets.isEmpty()) {
				throw new IOException("No target PSMs to combine");
			}

			File tmp = new File(combinedDir, "psms.tsv.tmp");
			OpenPsmTsvWriter writer = new OpenPsmTsvWriter(tmp, OpenPsmTsvWriter.pFindFormat);
			int written = 0;
			int noRt = 0;
			for (List<PFindPSM> list : List.of(targets, hosts)) {
				for (PFindPSM psm : list) {
					Double rt = rtMap.get(psm.getFileName() + "_" + psm.getScan());
					if (rt == null) {
						noRt++;
						continue;
					}
					psm.setRt(rt);
					writer.addPFindPSM(psm);
					written++;
				}
			}
			for (PFindPSM psm : decoys) {
				Double rt = rtMap.get(psm.getFileName() + "_" + psm.getScan());
				if (rt != null) {
					psm.setRt(rt);
				}
				writer.addPFindPSM(psm);
			}
			writer.close();
			Files.move(tmp.toPath(), combined.toPath(), StandardCopyOption.REPLACE_EXISTING);
			info("combined " + written + " target/host PSMs (" + hosts.size() + " host, " + noRt
					+ " dropped without a retention time) and " + decoys.size() + " decoys into "
					+ combined.getAbsolutePath());
			return combined;
		}

		// ------------------------------------------------------- 6. quant

		File quantify(File psms) throws IOException {
			File lfq = new File(combinedDir, "QuantifiedPeptides.tsv");
			if (hasContent(lfq)) {
				info("label-free quantification already exists: " + lfq.getAbsolutePath());
			} else {
				FlashLfqTask flashLfq = new FlashLfqTask(req.flashLfqExe.getAbsolutePath(), params.getThreadCount());
				boolean ok = flashLfq.run(params.getMetadata(), psms.getAbsolutePath(), rawDir.getAbsolutePath(), true,
						true);
				if (!ok || !hasContent(lfq)) {
					throw new IOException("FlashLFQ label-free quantification failed; expected " + lfq.getAbsolutePath());
				}
			}
			if (!MetaConstants.isobaricLabel.equals(quanMode)) {
				return lfq;
			}

			ArrayList<PFindPSM> all = new ArrayList<>();
			for (File sampleDir : sampleDirs) {
				File unfiltered = new File(new File(sampleDir, SAMPLE_DB_SEARCH_DIR), "pFind.spectra");
				if (unfiltered.isFile()) {
					all.addAll(Arrays.asList(new PFindPsmReader(unfiltered).getPsms()));
				}
			}
			PFindIsobaricTask isobaric = new PFindIsobaricTask(all.toArray(new PFindPSM[0]), spectraDir, combinedDir,
					params.getIsobaricTag(), params.getIsobaric(), params.getIsoCorFactor(), params.isCombineLabel(),
					params.getMetadata());
			Boolean ok = isobaric.runOnlyPeptide();
			File quanPep = isobaric.getQuanPepFile();
			if (ok == null || !ok || !hasContent(quanPep)) {
				throw new IOException("Isobaric peptide quantification failed; expected " + describe(quanPep));
			}
			info("isobaric quantification: " + quanPep.getAbsolutePath());
			return quanPep;
		}

		// ------------------------------------------------------- 8. variants

		void runPangenomeVariantCheck() {
			try {
				PangenomeVariantSettings settings = params.getPangenomeVariantSettings();
				File magVersionDir = mag.getCurrentFile();
				PangenomeVariantCheck.Request request = new PangenomeVariantCheck.Request();
				request.settings = settings;
				request.magVersionDir = magVersionDir;
				request.info = this::info;
				request.warn = msg -> info("WARNING: " + msg);
				if (!PangenomeVariantCheck.wouldRun(settings, magVersionDir)) {
					PangenomeVariantCheck.run(request);
					return;
				}
				ArrayList<File> spectraFiles = new ArrayList<>();
				for (File sampleDir : sampleDirs) {
					spectraFiles.add(new File(new File(sampleDir, SAMPLE_DB_SEARCH_DIR), FILTERED_SPECTRA));
				}
				double tolerance = params.getUndecidableToleranceDa();
				PFindVariantCandidateSource.Result source = PFindVariantCandidateSource.read(spectraFiles,
						mag.getIdentifier(), tolerance, this::info);
				info("pangenome variant check: pFind front end - " + source.describe());
				request.engine = VariantCandidate.Engine.PFIND;
				request.candidates = source.candidates;
				request.identifiedBasePeptides = source.basePeptides;
				request.searchedFastas = PFindVariantCandidateSource.searchedDatabases(sampleDirs,
						SAMPLE_DB_SEARCH_DIR, MARKER_SEARCH_DIR);
				request.decoyPrefix = PFindResultReader.REV;
				request.taxonomyDir = new File(resultRoot, "taxonomy_analysis");
				request.eggnogDir = mag.getOriginalFuncFolder();
				request.outputDir = combinedDir;
				request.undecidableToleranceDa = tolerance;
				request.fixedCysteineModDa = PFindVariantCandidateSource.FIXED_CYSTEINE_MOD_DA;
				PangenomeVariantCheck.run(request);
			} catch (RuntimeException e) {
				info("WARNING: pangenome variant check failed; workflow results are unaffected: " + e);
			}
		}

		// ------------------------------------------------------- helpers

		PFindTask newSearchTask(int pool) {
			PFindTask task = new PFindTask(pFindBin);
			task.setTotalThread(params.getThreadCount() / pool + 1);
			return task;
		}

		void addSearch(PFindTask task, int i, File outDir, File fasta) throws IOException {
			File pf2 = new File(pf2Files[i]);
			if (!pf2.isFile()) {
				throw new IOException("Spectra file not found: " + pf2.getAbsolutePath());
			}
			MetaOpenPar mop = params.getMop();
			task.searchPF2(new String[] { pf2.getAbsolutePath() }, outDir.getAbsolutePath(), fasta.getAbsolutePath(),
					params.getMs2ScanMode(), params.isOpenSearch(), true, mop.getOpenPsmFDR(),
					params.getFixedModifications().toArray(new String[0]),
					params.getVariableModifications().toArray(new String[0]), params.getEnzyme(),
					params.getMissCleavages(), params.getDigestMode(), params.getIsobaricTag());
		}

		File sampleDatabase(File sampleDir) {
			return new File(new File(sampleDir, MARKER_SEARCH_DIR), sampleDir.getName() + ".fasta");
		}

		int indexOf(String raw) {
			for (int i = 0; i < rawFiles.length; i++) {
				if (rawFiles[i].equals(raw)) {
					return i;
				}
			}
			throw new IllegalStateException(raw);
		}

		void step(String message) {
			info(message);
		}

		void info(String message) {
			cb.onMessage(TASK + ": " + message);
		}
	}

	/** pParse's .pf2 suffix for the MS2 scan mode, as the GUI derives it. */
	static String ms2Suffix(String ms2Type) {
		if (MetaConstants.HCD_FTMS.equals(ms2Type)) {
			return "_HCDFT";
		}
		if (MetaConstants.HCD_ITMS.equals(ms2Type)) {
			return "_HCDIT";
		}
		if (MetaConstants.CID_FTMS.equals(ms2Type)) {
			return "_CIDFT";
		}
		if (MetaConstants.CID_ITMS.equals(ms2Type)) {
			return "_CIDIT";
		}
		throw new IllegalArgumentException("Unknown MS2 scan mode '" + ms2Type + "'; expected one of "
				+ MetaConstants.HCD_FTMS + ", " + MetaConstants.HCD_ITMS + ", " + MetaConstants.CID_FTMS + ", "
				+ MetaConstants.CID_ITMS);
	}

	/** One search per 5 GB of physical memory, at least one - the GUI's pFind setting. */
	static int memoryBasedTaskCount() {
		long gb = ((com.sun.management.OperatingSystemMXBean) ManagementFactory.getOperatingSystemMXBean())
				.getTotalMemorySize() / 1024L / 1024L / 1024L;
		return (int) Math.max(1L, gb / 5L);
	}

	private static void copyLines(File in, PrintWriter out) throws IOException {
		try (BufferedReader reader = new BufferedReader(new FileReader(in))) {
			String line;
			while ((line = reader.readLine()) != null) {
				out.println(line);
			}
		}
	}

	private static boolean hasContent(File f) {
		return f != null && f.isFile() && f.length() > 0;
	}

	private static String describe(File f) {
		return f == null ? "null" : f.getAbsolutePath();
	}
}
