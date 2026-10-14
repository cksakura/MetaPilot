package uk.ac.quadram.metapilot.core.clean;

import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.LinkOption;
import java.nio.file.Path;
import java.nio.file.StandardOpenOption;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.Date;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.Set;

import org.json.JSONObject;

import uk.ac.quadram.metapilot.workflow.WorkflowType;

/**
 * Preview-first cleaner for MetaPilot result folders.
 */
public class ResultFolderCleaner {
	private static final boolean WORKFLOW_PROFILES_VERIFIED = true;

	public enum CleanMode {
		BASIC("Basic clean",
				"Keep final deliverables and workflow-critical files needed for rerun or resume, including parameters, logs, genome selection tables, manifests and refined search-space files."),
		DEEP("Deep clean",
				"Keep only final deliverables: root peptide/protein/summary files, taxonomy outputs, functional outputs, reports, and DIA combined outputs.");

		private final String displayName;
		private final String description;

		CleanMode(String displayName, String description) {
			this.displayName = displayName;
			this.description = description;
		}

		public String getDisplayName() {
			return displayName;
		}

		public String getDescription() {
			return description;
		}
	}

	public static final class DetectionResult {
		private final WorkflowType workflowType;
		private final String refinementMode;
		private final boolean highConfidence;
		private final List<String> evidence;
		private final List<String> warnings;

		public DetectionResult(WorkflowType workflowType, String refinementMode, boolean highConfidence,
				List<String> evidence, List<String> warnings) {
			this.workflowType = workflowType;
			this.refinementMode = refinementMode;
			this.highConfidence = highConfidence;
			this.evidence = Collections.unmodifiableList(new ArrayList<>(evidence));
			this.warnings = Collections.unmodifiableList(new ArrayList<>(warnings));
		}

		public WorkflowType getWorkflowType() {
			return workflowType;
		}

		public String getRefinementMode() {
			return refinementMode;
		}

		public boolean isHighConfidence() {
			return highConfidence;
		}

		public List<String> getEvidence() {
			return evidence;
		}

		public List<String> getWarnings() {
			return warnings;
		}

		public boolean hasWorkflow() {
			return workflowType != null;
		}
	}

	public static final class PreviewEntry {
		private final Path path;
		private final boolean directory;
		private final boolean keep;
		private final long sizeBytes;
		private final String reason;

		public PreviewEntry(Path path, boolean directory, boolean keep, long sizeBytes, String reason) {
			this.path = path;
			this.directory = directory;
			this.keep = keep;
			this.sizeBytes = sizeBytes;
			this.reason = reason;
		}

		public Path getPath() {
			return path;
		}

		public boolean isDirectory() {
			return directory;
		}

		public boolean isKeep() {
			return keep;
		}

		public long getSizeBytes() {
			return sizeBytes;
		}

		public String getReason() {
			return reason;
		}
	}

	public static final class CleanPreview {
		private final Path resultRoot;
		private final WorkflowType workflowType;
		private final String refinementMode;
		private final CleanMode mode;
		private final boolean removeConvertedSpectra;
		private final List<PreviewEntry> entries;
		private final List<String> warnings;
		private final long bytesToDelete;
		private final int filesToDelete;
		private final int directoriesToDelete;

		public CleanPreview(Path resultRoot, WorkflowType workflowType, String refinementMode, CleanMode mode,
				boolean removeConvertedSpectra, List<PreviewEntry> entries, List<String> warnings) {
			this.resultRoot = resultRoot;
			this.workflowType = workflowType;
			this.refinementMode = refinementMode;
			this.mode = mode;
			this.removeConvertedSpectra = removeConvertedSpectra;
			this.entries = Collections.unmodifiableList(new ArrayList<>(entries));
			this.warnings = Collections.unmodifiableList(new ArrayList<>(warnings));

			long bytes = 0L;
			int files = 0;
			int dirs = 0;
			for (PreviewEntry entry : entries) {
				if (!entry.isKeep()) {
					bytes += entry.getSizeBytes();
					if (entry.isDirectory()) {
						dirs++;
					} else {
						files++;
					}
				}
			}
			this.bytesToDelete = bytes;
			this.filesToDelete = files;
			this.directoriesToDelete = dirs;
		}

		public Path getResultRoot() {
			return resultRoot;
		}

		public WorkflowType getWorkflowType() {
			return workflowType;
		}

		public String getRefinementMode() {
			return refinementMode;
		}

		public CleanMode getMode() {
			return mode;
		}

		public boolean isRemoveConvertedSpectra() {
			return removeConvertedSpectra;
		}

		public List<PreviewEntry> getEntries() {
			return entries;
		}

		public List<String> getWarnings() {
			return warnings;
		}

		public long getBytesToDelete() {
			return bytesToDelete;
		}

		public int getFilesToDelete() {
			return filesToDelete;
		}

		public int getDirectoriesToDelete() {
			return directoriesToDelete;
		}
	}

	public static final class CleanResult {
		private final int deletedFiles;
		private final int deletedDirectories;
		private final long deletedBytes;
		private final Path logFile;

		public CleanResult(int deletedFiles, int deletedDirectories, long deletedBytes, Path logFile) {
			this.deletedFiles = deletedFiles;
			this.deletedDirectories = deletedDirectories;
			this.deletedBytes = deletedBytes;
			this.logFile = logFile;
		}

		public int getDeletedFiles() {
			return deletedFiles;
		}

		public int getDeletedDirectories() {
			return deletedDirectories;
		}

		public long getDeletedBytes() {
			return deletedBytes;
		}

		public Path getLogFile() {
			return logFile;
		}
	}

	private static final Set<String> FINAL_DIRECTORIES = Set.of("taxonomy", "taxonomy_analysis", "function",
			"functional_annotation", "report", "reports");
	private static final Set<String> WORKFLOW_RESULT_DIRECTORIES = Set.of("maxquant_search", "pfind_igc_search",
			"fragpipe_igc_search", "sage_igc_search", "pfind_mag_search", "fragpipe_mag_search", "sage_mag_search",
			"opensearch_mag_search", "marker_focus", "marker_comprehensive", "marker_deep", "seed_direct",
			"seed_focus", "seed_comprehensive", "seed_deep");
	private static final Set<String> ROOT_FINAL_FILES = Set.of("final_peptides.tsv", "final_proteins.tsv",
			"razor_proteins.tsv", "all_proteins.tsv", "final_summary.tsv", "summary.tsv", "proteingroups.txt",
			"peptides.txt", "summary.txt");
	private static final Set<String> DIA_COMBINED_DELIVERABLES = Set.of("combined.parquet", "combined.pr_matrix.tsv",
			"combined.pg_matrix.tsv", "combined.fasta", "combined.speclib", "combined.predicted.speclib",
			"combined.lib.parquet", "combined_report.tsv");
	private static final Set<String> TARGET_LIBRARY_EXTENSIONS = Set.of(".speclib", ".parquet", ".fasta", ".fas",
			".fa", ".faa", ".info", ".json", ".tsv");
	private static final Set<String> TEMP_FILE_SUFFIXES = Set.of(".pepindex", ".prot.xml", ".pepxml", ".tmp", ".temp",
			".bak");
	private static final Set<String> CONVERTED_SPECTRA_EXTENSIONS = Set.of(".mgf", ".mzml", ".mzxml", ".ms1", ".ms2",
			".pf1", ".pf2", ".dta");
	private static final Set<String> CONVERTED_SPECTRA_DIRECTORIES = Set.of("mgf", "mzml", "mzxml", "spectra",
			"converted_spectra", "converted spectra", "pfind_spectra", "pfind spectra");

	private static final Map<String, WorkflowType> PARAMETER_FILE_WORKFLOWS = createParameterFileWorkflowMap();

	public DetectionResult detect(Path resultRoot) {
		List<String> evidence = new ArrayList<>();
		List<String> warnings = new ArrayList<>();

		if (resultRoot == null || !Files.isDirectory(resultRoot)) {
			warnings.add("Select an existing result folder.");
			return new DetectionResult(null, null, false, evidence, warnings);
		}

		Map<WorkflowType, Integer> scores = new LinkedHashMap<>();
		String refinementMode = null;
		String rootName = resultRoot.getFileName() == null ? "" : resultRoot.getFileName().toString().toLowerCase(Locale.ROOT);
		WorkflowType fromRootName = workflowFromRootName(rootName);
		if (fromRootName != null) {
			addScore(scores, fromRootName, 5);
			evidence.add("Result folder name matches workflow root: " + rootName);
		}

		try {
			List<Path> candidates = new ArrayList<>();
			try (var stream = Files.walk(resultRoot, 2)) {
				stream.filter(path -> Files.isRegularFile(path, LinkOption.NOFOLLOW_LINKS)).forEach(candidates::add);
			}

			for (Path candidate : candidates) {
				String name = candidate.getFileName().toString().toLowerCase(Locale.ROOT);
				WorkflowType fromName = PARAMETER_FILE_WORKFLOWS.get(name);
				if (fromName != null) {
					addScore(scores, fromName, 2);
					evidence.add("Found parameter file: " + resultRoot.relativize(candidate));
				}

				if (name.endsWith(".json") && Files.size(candidate) < 2_000_000L) {
					JsonDetection jsonDetection = detectFromJson(candidate);
					if (jsonDetection.workflowType != null) {
						addScore(scores, jsonDetection.workflowType, 5);
						evidence.add("Read workflow '" + jsonDetection.workflowType.getAnalysisType() + "' from "
								+ resultRoot.relativize(candidate));
					}
					if (refinementMode == null && jsonDetection.refinementMode != null) {
						refinementMode = jsonDetection.refinementMode;
					}
				}
			}
		} catch (Exception e) {
			warnings.add("Workflow auto-detection could not read every candidate file: " + e.getMessage());
		}

		WorkflowType workflowType = null;
		int bestScore = 0;
		boolean tie = false;
		for (Map.Entry<WorkflowType, Integer> entry : scores.entrySet()) {
			if (entry.getValue() > bestScore) {
				workflowType = entry.getKey();
				bestScore = entry.getValue();
				tie = false;
			} else if (entry.getValue() == bestScore) {
				tie = true;
			}
		}

		if (workflowType == null) {
			warnings.add("Workflow could not be detected automatically. Select it manually before previewing.");
		} else if (tie) {
			warnings.add("Multiple workflow markers were found. Review the selected workflow before cleaning.");
		}

		return new DetectionResult(workflowType, refinementMode, workflowType != null && !tie && bestScore >= 5,
				evidence, warnings);
	}

	public CleanPreview preview(Path selectedRoot, WorkflowType workflowType, String refinementMode, CleanMode mode)
			throws IOException {
		return preview(selectedRoot, workflowType, refinementMode, mode, false);
	}

	public CleanPreview preview(Path selectedRoot, WorkflowType workflowType, String refinementMode, CleanMode mode,
			boolean removeConvertedSpectra) throws IOException {
		if (selectedRoot == null) {
			throw new IOException("Result folder is not selected.");
		}
		if (mode == null) {
			throw new IOException("Clean mode is not selected.");
		}
		if (workflowType == null) {
			throw new IOException("Workflow is not selected. Select the exact workflow before previewing.");
		}
		refuseUntilWorkflowProfilesVerified();
		Path resultRoot = selectedRoot.toRealPath();
		if (!Files.isDirectory(resultRoot)) {
			throw new IOException("Result folder does not exist: " + resultRoot);
		}

		List<String> warnings = new ArrayList<>();
		if (removeConvertedSpectra) {
			warnings.add(
					"Converted spectra are selected for removal. Reruns will need to perform spectra conversion again.");
		}

		List<Path> allPaths = new ArrayList<>();
		try (var stream = Files.walk(resultRoot)) {
			stream.filter(path -> !path.equals(resultRoot)).forEach(allPaths::add);
		}
		allPaths.sort(Comparator.comparing(path -> resultRoot.relativize(path).toString()));
		validateSelectedRoot(resultRoot, workflowType, refinementMode, allPaths);

		Set<Path> keptFiles = new HashSet<>();
		Map<Path, String> reasons = new LinkedHashMap<>();
		for (Path path : allPaths) {
			if (Files.isDirectory(path, LinkOption.NOFOLLOW_LINKS)) {
				continue;
			}
			KeepDecision decision = shouldKeepFile(resultRoot, path, workflowType, mode, removeConvertedSpectra);
			if (decision.keep) {
				keptFiles.add(path);
				reasons.put(path, decision.reason);
			}
		}

		List<PreviewEntry> entries = new ArrayList<>();
		for (Path path : allPaths) {
			boolean directory = Files.isDirectory(path, LinkOption.NOFOLLOW_LINKS);
			boolean keep;
			String reason;
			if (directory) {
				KeepDecision directoryDecision = shouldKeepDirectory(resultRoot, path, workflowType, mode);
				boolean containsKeptFile = keptFiles.stream().anyMatch(file -> file.startsWith(path));
				keep = directoryDecision.keep || containsKeptFile;
				reason = directoryDecision.keep ? directoryDecision.reason
						: containsKeptFile ? "Contains retained files" : "Not retained by selected clean mode";
			} else {
				keep = keptFiles.contains(path);
				reason = reasons.getOrDefault(path, "Not retained by selected clean mode");
			}
			long size = directory ? 0L : safeSize(path);
			entries.add(new PreviewEntry(path, directory, keep, size, reason));
		}

		return new CleanPreview(resultRoot, workflowType, refinementMode, mode, removeConvertedSpectra, entries,
				warnings);
	}

	public CleanResult clean(CleanPreview preview) throws IOException {
		refuseUntilWorkflowProfilesVerified();
		if (preview == null) {
			throw new IOException("No clean preview is available.");
		}

		Path root = preview.getResultRoot().toRealPath();
		List<PreviewEntry> deleteEntries = new ArrayList<>();
		for (PreviewEntry entry : preview.getEntries()) {
			if (!entry.isKeep()) {
				deleteEntries.add(entry);
			}
		}
		deleteEntries.sort((a, b) -> b.getPath().getNameCount() - a.getPath().getNameCount());

		int files = 0;
		int dirs = 0;
		long bytes = 0L;
		List<String> logLines = new ArrayList<>();
		logLines.add("MetaPilot result folder clean");
		logLines.add("Generated: " + new SimpleDateFormat("yyyy-MM-dd HH:mm:ss").format(new Date()));
		logLines.add("Root: " + root);
		logLines.add("Workflow: "
				+ (preview.getWorkflowType() != null ? preview.getWorkflowType().getAnalysisType() : "UNKNOWN"));
		logLines.add(
				"Refinement mode: " + (preview.getRefinementMode() != null ? preview.getRefinementMode() : "UNKNOWN"));
		logLines.add("Mode: " + preview.getMode().getDisplayName());
		logLines.add("Remove converted spectra: " + preview.isRemoveConvertedSpectra());
		if (!preview.getWarnings().isEmpty()) {
			logLines.add("");
			logLines.add("Warnings:");
			for (String warning : preview.getWarnings()) {
				logLines.add("- " + warning);
			}
		}
		logLines.add("");
		logLines.add("Deleted paths:");

		for (PreviewEntry entry : deleteEntries) {
			Path path = entry.getPath().toAbsolutePath().normalize();
			if (!path.startsWith(root) || path.equals(root)) {
				throw new IOException("Refusing to delete outside the selected result folder: " + path);
			}
			if (!Files.exists(path, LinkOption.NOFOLLOW_LINKS)) {
				continue;
			}
			if (entry.isDirectory()) {
				try {
					Files.delete(path);
					dirs++;
					logLines.add("[DIR] " + root.relativize(path));
				} catch (IOException e) {
					try (var children = Files.list(path)) {
						if (children.findAny().isPresent()) {
							continue;
						}
					}
					throw e;
				}
			} else {
				long size = entry.getSizeBytes();
				Files.delete(path);
				files++;
				bytes += size;
				logLines.add("[FILE] " + root.relativize(path) + " (" + size + " bytes)");
			}
		}

		Path logFile = root.resolve("metapilot_clean_log.txt");
		logLines.add("");
		logLines.add("Deleted files: " + files);
		logLines.add("Deleted directories: " + dirs);
		logLines.add("Deleted bytes: " + bytes);
		logLines.add("");
		logLines.add("Retained paths:");
		for (PreviewEntry entry : preview.getEntries()) {
			if (entry.isKeep()) {
				String type = entry.isDirectory() ? "DIR" : "FILE";
				logLines.add("[" + type + "] " + root.relativize(entry.getPath()) + " - " + entry.getReason());
			}
		}
		Files.write(logFile, logLines, StandardCharsets.UTF_8, StandardOpenOption.CREATE,
				StandardOpenOption.TRUNCATE_EXISTING);

		return new CleanResult(files, dirs, bytes, logFile);
	}

	private KeepDecision shouldKeepDirectory(Path root, Path path, WorkflowType workflowType, CleanMode mode) {
		Path relative = root.relativize(path);
		if (Files.isSymbolicLink(path)) {
			return KeepDecision.keep("Symbolic link retained; cleaner never follows or deletes symlinks");
		}
		if (containsSegment(relative, FINAL_DIRECTORIES)) {
			return KeepDecision.keep("Final taxonomy/function/report folder");
		}
		if (mode == CleanMode.BASIC) {
			return KeepDecision.keep("Basic clean retains workflow state directories unless all children are known disposable files");
		}
		if (isDiaCombinedDirectory(workflowType, relative)) {
			return KeepDecision.keep("DIA combined deliverable folder");
		}
		return KeepDecision.delete();
	}

	private KeepDecision shouldKeepFile(Path root, Path path, WorkflowType workflowType, CleanMode mode,
			boolean removeConvertedSpectra) {
		Path relative = root.relativize(path);
		String name = path.getFileName().toString().toLowerCase(Locale.ROOT);

		if (Files.isSymbolicLink(path)) {
			return KeepDecision.keep("Symbolic link retained; cleaner never follows or deletes symlinks");
		}
		if (name.equals("metapilot_clean_log.txt")) {
			return KeepDecision.keep("Clean history");
		}
		if (isConvertedSpectraFile(relative, name)) {
			if (removeConvertedSpectra) {
				return KeepDecision.delete("Converted spectra selected for removal; rerun will reconvert spectra");
			}
			return KeepDecision.keep("Converted spectra retained for rerun/resume");
		}
		if (containsSegment(relative, FINAL_DIRECTORIES)) {
			return KeepDecision.keep("Final taxonomy/function/report output");
		}
		if (isFinalDeliverableFile(workflowType, relative, name)) {
			return KeepDecision.keep("Top-level peptide/protein/summary/report output");
		}
		if (isDiaCombinedDeliverable(workflowType, relative, name)) {
			return KeepDecision.keep("DIA combined matrix/library deliverable");
		}
		boolean pfindMagDeepClean = workflowType == WorkflowType.DDA_MAG_PFIND && mode == CleanMode.DEEP;
		boolean fragpipeMagDeepClean = workflowType == WorkflowType.DDA_MAG_FRAGPIPE && mode == CleanMode.DEEP;
		if (!pfindMagDeepClean && !fragpipeMagDeepClean && isDdaCombinedDeliverable(workflowType, relative, name)) {
			return KeepDecision.keep("DDA combined search result table");
		}
		if (pfindMagDeepClean && isPfindMagDeepCombinedFile(relative, name)) {
			return KeepDecision.keep("pFind combined/ resume-gate file retained in deep clean "
					+ "(combineResults()/lfQuant() check this file's existence to skip re-combining/re-quantifying on rerun; "
					+ "other combined/ outputs are regenerated from psms.tsv, so deep clean discards them)");
		}
		if (fragpipeMagDeepClean && isFragpipeMagDeepCombinedFile(relative, name)) {
			return KeepDecision.keep("FragPipe combined/ resume-gate file retained in deep clean "
					+ "(runInternal() checks combined_modified_peptide.tsv's existence to skip the combined FragPipe search, "
					+ "and combined.fasta's existence to skip rebuilding it from every sample's protein.fas; "
					+ "combined_protein.tsv is never read back by MetaPilot, so deep clean discards it)");
		}
		if (isTargetTopLevelAsset(workflowType, relative, name)) {
			return KeepDecision.keep("DIA target library/input asset retained");
		}
		if (isPfindMagPerSampleDeepDeliverable(workflowType, relative, name)) {
			return KeepDecision.keep(
					"pFind per-sample search result retained (used to build combined results / rerun; kept even in deep clean)");
		}
		if (isFragpipeMagPerSampleDeepDeliverable(workflowType, relative, name)) {
			return KeepDecision.keep(
					"FragPipe per-sample search result retained (used to build combined results / rerun; kept even in deep clean)");
		}

		// Open-search workflow: without this the DEEP fallback below
		// (KeepDecision.delete()) would destroy its resume sentinels and the
		// validation intermediates, which no generic rule recognises.
		if (workflowType == WorkflowType.DDA_MAG_OPENSEARCH
				&& uk.ac.quadram.metapilot.core.opensearch.OpenSearchLayout.deepCleanKeepFileNames().contains(name)) {
			return KeepDecision.keep("Open-search resume gate / validation output retained even in deep clean "
					+ "(the runner skips a step when its sentinel exists, so deleting one of these would leave the "
					+ "step skipped and its output never regenerated)");
		}
		if (workflowType == WorkflowType.DDA_MAG_OPENSEARCH && name.endsWith(".pepXML")) {
			return KeepDecision.keep("Open-search per-sample pepXML retained (stage-1/stage-2 search resume gates)");
		}
		// Crystal-C's own artifacts, kept for the same provenance reason
		// fragger_open.params is: the refined _c.pepXML beside them is a resume
		// sentinel, so deleting the report and parameters leaves a refinement that
		// will be skipped on the next run with no record of what it did.
		if (workflowType == WorkflowType.DDA_MAG_OPENSEARCH
				&& (name.endsWith(uk.ac.quadram.metapilot.core.crystalc.CrystalCHeadlessRunner.REPORT_SUFFIX)
						|| (name.startsWith("crystalc-") && name.endsWith(".params")))) {
			return KeepDecision.keep("Crystal-C precursor-refinement report / parameters retained even in deep "
					+ "clean (provenance for the _c.pepXML the validation reads)");
		}

		if (mode == CleanMode.BASIC) {
			if (isOpenSearchMagGenomeSelectionAudit(workflowType, relative, name)) {
				return KeepDecision.keep("Open-search genome-selection audit artifact retained for basic clean");
			}
			if (isPfindMagGenomeSelectionAudit(workflowType, relative, name)) {
				return KeepDecision.keep("pFind genome-selection audit artifact retained for basic clean");
			}
			if (isFragpipeMagGenomeSelectionAudit(workflowType, relative, name)) {
				return KeepDecision.keep("FragPipe genome-selection audit artifact retained for basic clean");
			}
			if (isExplicitTemporaryFile(relative, name)) {
				return KeepDecision.delete("Audited temporary/cache file not required for rerun");
			}
			return KeepDecision.keep("Basic clean retains unclassified workflow state to preserve rerun/resume safety");
		}

		return KeepDecision.delete();
	}

	private static boolean isFinalDeliverableFile(WorkflowType workflowType, Path relative, String name) {
		// QPX exports ({PREFIX}.psm.parquet, written beside final_peptides.tsv) are final
		// deliverables too.
		if (!ROOT_FINAL_FILES.contains(name) && !name.startsWith("final_") && !name.endsWith(".psm.parquet")) {
			return false;
		}
		if (relative.getNameCount() == 1) {
			return true;
		}
		if (relative.getNameCount() == 2) {
			String parent = relative.getName(0).toString().toLowerCase(Locale.ROOT);
			return WORKFLOW_RESULT_DIRECTORIES.contains(parent)
					|| (workflowType == WorkflowType.DIA_REALTIME && parent.equals("combined"))
					|| (workflowType == WorkflowType.DDA_MAG_OPENSEARCH && uk.ac.quadram.metapilot.core.opensearch.OpenSearchLayout
							.modeFolderNames().contains(parent));
		}
		return false;
	}

	private static boolean isDiaCombinedDirectory(WorkflowType workflowType, Path relative) {
		return isDiaWorkflow(workflowType) && relative.getNameCount() == 1
				&& firstSegment(relative).equals("combined");
	}

	private static boolean isDiaCombinedDeliverable(WorkflowType workflowType, Path relative, String name) {
		if (!isDiaWorkflow(workflowType)) {
			return false;
		}
		if (relative.getNameCount() == 2 && firstSegment(relative).equals("combined")) {
			return DIA_COMBINED_DELIVERABLES.contains(name) || ROOT_FINAL_FILES.contains(name)
					|| name.startsWith("final_");
		}
		if (workflowType == WorkflowType.DIA_SEED && name.matches("round\\d+\\.(pr|pg)_matrix\\.tsv")) {
			return true;
		}
		return false;
	}

	private static boolean isDdaCombinedDeliverable(WorkflowType workflowType, Path relative, String name) {
		if (workflowType == null || isDiaWorkflow(workflowType) || relative.getNameCount() < 2
				|| !firstSegment(relative).equals("combined")) {
			return false;
		}
		return name.endsWith(".tsv") || name.endsWith(".txt");
	}

	/**
	 * The only two files under a pFind DDA MAG {@code combined/} folder that
	 * gate resume/rerun: {@code combineResults()} skips re-merging every
	 * sample's PSMs when {@code combined/psms.tsv} already exists (and is
	 * over 200 bytes), and {@code lfQuant()} skips re-running FlashLFQ
	 * label-free quantification when {@code combined/QuantifiedPeptides.tsv}
	 * already exists (and is non-empty). Every other file under
	 * {@code combined/} (`id.tsv`, `QuantifiedProteins.tsv`,
	 * `QuantifiedPeaks.tsv`, `ms2_intensity.tsv`, the isobaric
	 * `IsobaricQuantified*.tsv` outputs) is regenerated from
	 * {@code psms.tsv} on the next run and is not itself checked by any
	 * resume logic, so deep clean discards them - only these two survive.
	 * Basic clean is unaffected: it still keeps every {@code combined/}
	 * file via the generic {@link #isDdaCombinedDeliverable} rule.
	 */
	private static boolean isPfindMagDeepCombinedFile(Path relative, String name) {
		if (relative.getNameCount() < 2 || !firstSegment(relative).equals("combined")) {
			return false;
		}
		return name.equals("psms.tsv") || name.equals("quantifiedpeptides.tsv");
	}

	private static boolean isConvertedSpectraFile(Path relative, String name) {
		if (name.endsWith(".raw.dia") || name.endsWith(".wiff.dia")) {
			return true;
		}
		boolean inConvertedSpectraDirectory = false;
		for (int i = 0; i < relative.getNameCount() - 1; i++) {
			String segment = relative.getName(i).toString().toLowerCase(Locale.ROOT);
			if (CONVERTED_SPECTRA_DIRECTORIES.contains(segment)) {
				inConvertedSpectraDirectory = true;
				break;
			}
		}
		if (inConvertedSpectraDirectory) {
			return true;
		}
		for (String extension : CONVERTED_SPECTRA_EXTENSIONS) {
			if (!extension.equals(".mzml") && !extension.equals(".mzxml") && name.endsWith(extension)) {
				return true;
			}
		}
		return false;
	}

	private static boolean isTargetTopLevelAsset(WorkflowType workflowType, Path relative, String name) {
		if (workflowType != WorkflowType.DIA_TARGET || relative.getNameCount() != 1) {
			return false;
		}
		for (String extension : TARGET_LIBRARY_EXTENSIONS) {
			if (name.endsWith(extension)) {
				return true;
			}
		}
		return false;
	}

	/**
	 * pFind DDA MAG per-sample search results that later pipeline steps read
	 * back within the same run (marker-DB search feeds sample-specific FASTA
	 * extraction; sample-DB search feeds PSM combination, isobaric
	 * quantification, and the {@code exportReport()} sample-list fallback) -
	 * see {@code DDAMagPfindTask}. These are kept even in deep clean, unlike
	 * other per-sample workflow state, because deleting them would silently
	 * force a full re-search on any future rerun/resume of this result
	 * folder. Matched by the last one or two path segments so both the
	 * current {@code search_results/<sample>/...} layout and any
	 * not-yet-migrated legacy {@code <sample>/...} layout are recognized -
	 * see {@code MagResultLayoutUtil}.
	 */
	private static boolean isPfindMagPerSampleDeepDeliverable(WorkflowType workflowType, Path relative, String name) {
		if (workflowType != WorkflowType.DDA_MAG_PFIND) {
			return false;
		}
		int count = relative.getNameCount();
		if (count < 2) {
			return false;
		}
		String parent = relative.getName(count - 2).toString().toLowerCase(Locale.ROOT);
		if (parent.equals("marker")) {
			if (name.equals("pfind-filtered.spectra")) {
				return true;
			}
			if (count >= 3 && name.endsWith(".fasta")) {
				String sampleFolder = relative.getName(count - 3).toString().toLowerCase(Locale.ROOT);
				return name.equals(sampleFolder + ".fasta");
			}
			return false;
		}
		if (parent.equals("sample_db")) {
			return name.equals("pfind-filtered.spectra") || name.equals("pfind.spectra");
		}
		if (parent.equals("dbreducer") && count >= 3) {
			String grandparent = relative.getName(count - 3).toString().toLowerCase(Locale.ROOT);
			if (grandparent.equals("sample_db")) {
				return name.equals("dbreducer.fasta");
			}
		}
		return false;
	}

	/**
	 * pFind DDA MAG genome-selection audit artifacts ({@code marker/*_evidence.tsv},
	 * {@code marker/*_ranked.tsv}, {@code marker/*_selected.tsv}). Nothing
	 * reads these back later, but they are diagnostic records of why each
	 * genome was selected into a sample's search database, so they are kept
	 * explicitly for basic clean (basic clean's generic "keep unclassified"
	 * fallback would already retain them today, but this makes the
	 * retention explicit and independent of that default). Deep clean still
	 * discards them - only {@link #isPfindMagPerSampleDeepDeliverable} files
	 * survive deep clean.
	 */
	private static boolean isPfindMagGenomeSelectionAudit(WorkflowType workflowType, Path relative, String name) {
		if (workflowType != WorkflowType.DDA_MAG_PFIND) {
			return false;
		}
		int count = relative.getNameCount();
		if (count < 2) {
			return false;
		}
		String parent = relative.getName(count - 2).toString().toLowerCase(Locale.ROOT);
		if (!parent.equals("marker")) {
			return false;
		}
		return name.endsWith("_evidence.tsv") || name.endsWith("_ranked.tsv") || name.endsWith("_selected.tsv");
	}

	/**
	 * FragPipe DDA MAG per-sample search results that later pipeline steps
	 * read back within the same run, or that gate resume of a sub-step - see
	 * {@code MagFragPipeHeadlessRunner}. {@code marker/<expFolder>_1/psm.tsv}
	 * gates {@code runStepAParallel()}'s per-sample resume-skip and is a
	 * required input to Step B's genome selection on every run.
	 * {@code marker/<sample>.fasta} gates Step B's FASTA-rebuild sub-step and
	 * is the search database Step C needs if it has to rerun.
	 * {@code sample_db/<expFolder>_1/protein.fas} gates {@code
	 * runStepCParallel()}'s per-sample resume-skip and is a required input
	 * to build {@code combined/combined.fasta}. Note {@code psm.tsv} and
	 * {@code protein.fas} sit one level deeper than pFind's equivalent files
	 * - inside FragPipe's own per-experiment output folder
	 * ({@code marker/<expFolder>_1/...}, not directly in {@code marker/...})
	 * - so this matches by grandparent segment for those two, and by parent
	 * segment (like pFind) for the FASTA. Deleting any of these would force
	 * a full re-search on any future rerun/resume of this result folder.
	 */
	private static boolean isFragpipeMagPerSampleDeepDeliverable(WorkflowType workflowType, Path relative,
			String name) {
		if (workflowType != WorkflowType.DDA_MAG_FRAGPIPE) {
			return false;
		}
		int count = relative.getNameCount();
		if (count < 2) {
			return false;
		}
		String parent = relative.getName(count - 2).toString().toLowerCase(Locale.ROOT);
		if (parent.equals("marker") && count >= 3 && name.endsWith(".fasta")) {
			String sampleFolder = relative.getName(count - 3).toString().toLowerCase(Locale.ROOT);
			if (name.equals(sampleFolder + ".fasta")) {
				return true;
			}
		}
		if (count >= 3) {
			String grandparent = relative.getName(count - 3).toString().toLowerCase(Locale.ROOT);
			if (grandparent.equals("marker") && name.equals("psm.tsv")) {
				return true;
			}
			if (grandparent.equals("sample_db") && name.equals("protein.fas")) {
				return true;
			}
		}
		return false;
	}

	/**
	 * FragPipe DDA MAG genome-selection audit artifacts
	 * ({@code marker/*_evidence.tsv}, {@code marker/*_ranked.tsv},
	 * {@code marker/*_selected.tsv}, written with a {@code marker_genome_}
	 * or {@code genome_} prefix depending on selector mode). Unlike pFind,
	 * {@code selectGenomes()} is not itself resume-gated - it recomputes
	 * these every time Step B runs for a sample, regardless of whether they
	 * already exist - but nothing reads them back afterward either, so they
	 * are still just diagnostic records. Kept explicitly for basic clean;
	 * deep clean discards them - only
	 * {@link #isFragpipeMagPerSampleDeepDeliverable} files survive deep clean.
	 */
	/**
	 * Genome-selection audit artifacts for the open-search workflow, using the
	 * same {@code *_evidence/_ranked/_selected.tsv} convention as the other MAG
	 * workflows so the shape of this rule carries over.
	 */
	private static boolean isOpenSearchMagGenomeSelectionAudit(WorkflowType workflowType, Path relative, String name) {
		if (workflowType != WorkflowType.DDA_MAG_OPENSEARCH) {
			return false;
		}
		for (String suffix : uk.ac.quadram.metapilot.core.opensearch.OpenSearchLayout
				.genomeSelectionAuditSuffixes()) {
			if (name.endsWith(suffix)) {
				return true;
			}
		}
		return false;
	}

	private static boolean isFragpipeMagGenomeSelectionAudit(WorkflowType workflowType, Path relative, String name) {
		if (workflowType != WorkflowType.DDA_MAG_FRAGPIPE) {
			return false;
		}
		int count = relative.getNameCount();
		if (count < 2) {
			return false;
		}
		String parent = relative.getName(count - 2).toString().toLowerCase(Locale.ROOT);
		if (!parent.equals("marker")) {
			return false;
		}
		return name.endsWith("_evidence.tsv") || name.endsWith("_ranked.tsv") || name.endsWith("_selected.tsv");
	}

	/**
	 * Files under a FragPipe DDA MAG {@code combined/} folder that survive
	 * deep clean: {@code runInternal()} skips the combined FragPipe search
	 * (the expensive step) when {@code combined/combined_modified_peptide.tsv}
	 * already exists and is non-empty, and skips rebuilding the combined
	 * FASTA - which requires re-reading every sample's {@code protein.fas} -
	 * when {@code combined/combined.fasta} already exists and is non-empty.
	 * {@code combined/<expName>/psm.tsv} - one such subfolder per experiment
	 * in the combined manifest, named after that experiment's intensity
	 * column title in {@code combined_modified_peptide.tsv} - is not itself
	 * resume-gated, but {@code FragpipePepReader.getPSMs()} (used by
	 * {@code FragpipePostAnalysisTsvExporter}, called every time Step G
	 * post-analysis runs, including when Step D/E/F's combined search was
	 * itself skipped via the {@code combined_modified_peptide.tsv} resume
	 * check above) reads PSM-level fields (Hyperscore, PeptideProphet
	 * Probability, Assigned Modifications) from these files -
	 * {@code combined_modified_peptide.tsv} is peptide-level only and has no
	 * PSM-level detail of its own. Deleting them produces a "peptides=N,
	 * psms=0" result on the next run without ever re-running the combined
	 * search to regenerate them, since that step's own resume check
	 * (checking {@code combined_modified_peptide.tsv}, not these files)
	 * still finds a reason to skip. {@code combined/combined_protein.tsv} is
	 * written by the same combined search but is never read back anywhere
	 * in MetaPilot (confirmed by grep across the codebase), so unlike
	 * pFind's {@code QuantifiedPeptides.tsv} it does not gate anything and
	 * deep clean discards it; the root-level {@code final_peptides.tsv}/
	 * {@code razor_proteins.tsv}/{@code all_proteins.tsv} deliverables are
	 * built from {@code combined_modified_peptide.tsv} plus these per-experiment
	 * {@code psm.tsv} files. {@code combined/raws.fp-manifest} and
	 * {@code combined/combined.workflow} are unconditionally rewritten only
	 * when the combined search actually runs (never on the skip path), so
	 * they are also safe to discard. Basic clean is unaffected: it still
	 * keeps every {@code combined/} file via the generic
	 * {@link #isDdaCombinedDeliverable} rule.
	 */
	private static boolean isFragpipeMagDeepCombinedFile(Path relative, String name) {
		if (relative.getNameCount() < 2 || !firstSegment(relative).equals("combined")) {
			return false;
		}
		if (relative.getNameCount() >= 3 && name.equals("psm.tsv")) {
			return true;
		}
		return name.equals("combined_modified_peptide.tsv") || name.equals("combined.fasta");
	}

	private static boolean isExplicitTemporaryFile(Path relative, String name) {
		if (name.equals("combined.prot.xml") || name.endsWith("_uncalibrated.mzml")) {
			return true;
		}
		for (String suffix : TEMP_FILE_SUFFIXES) {
			if (name.endsWith(suffix)) {
				return true;
			}
		}
		if (relative.getNameCount() > 1) {
			String parent = relative.getName(relative.getNameCount() - 2).toString().toLowerCase(Locale.ROOT);
			return parent.equals("temp") || parent.equals("tmp");
		}
		return false;
	}

	private static void validateSelectedRoot(Path resultRoot, WorkflowType workflowType, String refinementMode,
			List<Path> allPaths) throws IOException {
		if (workflowType == null) {
			throw new IOException("Workflow is not selected. Select the exact workflow before previewing.");
		}
		if (workflowType == WorkflowType.DIA_DEEP_COMPREHENSIVE) {
			throw new IOException("DIA Deep Comprehensive cleaning is not exposed in the main app and has no "
					+ "standalone cleaner profile. Refusing to preview.");
		}
		String rootName = resultRoot.getFileName() == null ? "" : resultRoot.getFileName().toString().toLowerCase(Locale.ROOT);
		List<String> nestedWorkflowRoots = nestedWorkflowRoots(resultRoot, allPaths);
		if (!nestedWorkflowRoots.isEmpty()) {
			throw new IOException("The selected folder looks like a parent output folder, not a single workflow result folder. "
					+ "Nested workflow result folder(s) detected: " + String.join(", ", nestedWorkflowRoots)
					+ ". Select one exact workflow result folder before previewing, for example: "
					+ nestedWorkflowRoots.get(0) + ".");
		}
		Set<String> allowedRoots = allowedRootDirectories(workflowType, refinementMode);
		if (!allowedRoots.isEmpty() && !allowedRoots.contains(rootName)) {
			throw new IOException("The selected folder name '" + resultRoot.getFileName()
					+ "' does not match the selected workflow "
					+ (workflowType != null ? workflowType.getAnalysisType() : "UNKNOWN")
					+ ". Expected one of: " + String.join(", ", allowedRoots)
					+ ". Select the exact workflow result folder or choose the correct workflow.");
		}
		if (allowedRoots.isEmpty() && WORKFLOW_RESULT_DIRECTORIES.contains(rootName)) {
			throw new IOException("The selected folder name '" + resultRoot.getFileName()
					+ "' looks like a workflow result folder, but it does not match the selected workflow "
					+ (workflowType != null ? workflowType.getAnalysisType() : "UNKNOWN")
					+ ". Choose the correct workflow before previewing.");
		}
		if (workflowType == WorkflowType.DIA_TARGET && !hasDiaTargetEvidence(resultRoot, allPaths)) {
			throw new IOException("The selected folder does not contain enough DIA Target evidence. Expected a target "
					+ "library asset at the selected root or a combined DIA output folder. Refusing to preview.");
		}
		if (workflowType == WorkflowType.DIA_REALTIME && !hasDiaRealtimeEvidence(resultRoot, allPaths)) {
			throw new IOException("The selected folder does not contain enough DIA Real-time evidence. Expected a "
					+ "combined output folder or per-file DIA parquet folders. Refusing to preview.");
		}
	}

	private static List<String> nestedWorkflowRoots(Path resultRoot, List<Path> allPaths) {
		String rootName = resultRoot.getFileName() == null ? "" : resultRoot.getFileName().toString().toLowerCase(Locale.ROOT);
		if (WORKFLOW_RESULT_DIRECTORIES.contains(rootName)) {
			return Collections.emptyList();
		}
		List<String> nestedWorkflowRoots = new ArrayList<>();
		for (Path path : allPaths) {
			Path relative = resultRoot.relativize(path);
			if (relative.getNameCount() == 1 && Files.isDirectory(path, LinkOption.NOFOLLOW_LINKS)) {
				String name = relative.getFileName().toString().toLowerCase(Locale.ROOT);
				if (WORKFLOW_RESULT_DIRECTORIES.contains(name)) {
					nestedWorkflowRoots.add(relative.toString());
				}
			}
		}
		return nestedWorkflowRoots;
	}

	private static Set<String> allowedRootDirectories(WorkflowType workflowType, String refinementMode) {
		if (workflowType == null) {
			return Collections.emptySet();
		}
		switch (workflowType) {
		case DDA_IGC_MAXQUANT:
			return Set.of("maxquant_search");
		case DDA_IGC_PFIND:
			return Set.of("pfind_igc_search");
		case DDA_IGC_FRAGPIPE:
			return Set.of("fragpipe_igc_search");
		case DDA_IGC_SAGE:
			return Set.of("sage_igc_search");
		case DDA_MAG_PFIND:
			return Set.of("pfind_mag_search");
		case DDA_MAG_FRAGPIPE:
			return Set.of("fragpipe_mag_search");
		case DDA_MAG_SAGE:
			return Set.of("sage_mag_search");
		case DDA_MAG_OPENSEARCH:
			return Set.of(uk.ac.quadram.metapilot.core.opensearch.OpenSearchLayout.ENGINE_SEARCH_DIR_KEY);
		case DIA_MARKER:
			return allowedDiaMarkerRoots(refinementMode);
		case DIA_SEED:
			return allowedDiaSeedRoots(refinementMode);
		default:
			return Collections.emptySet();
		}
	}

	private static boolean hasDiaTargetEvidence(Path resultRoot, List<Path> allPaths) {
		for (Path path : allPaths) {
			Path relative = resultRoot.relativize(path);
			String name = path.getFileName() == null ? "" : path.getFileName().toString().toLowerCase(Locale.ROOT);
			if (relative.getNameCount() == 1 && Files.isRegularFile(path, LinkOption.NOFOLLOW_LINKS)
					&& isTargetTopLevelAsset(WorkflowType.DIA_TARGET, relative, name)) {
				return true;
			}
			if (relative.getNameCount() == 1 && Files.isDirectory(path, LinkOption.NOFOLLOW_LINKS)
					&& firstSegment(relative).equals("combined")) {
				return true;
			}
			if (relative.getNameCount() == 2 && firstSegment(relative).equals("combined")
					&& DIA_COMBINED_DELIVERABLES.contains(name)) {
				return true;
			}
		}
		return false;
	}

	private static boolean hasDiaRealtimeEvidence(Path resultRoot, List<Path> allPaths) {
		for (Path path : allPaths) {
			Path relative = resultRoot.relativize(path);
			String name = path.getFileName() == null ? "" : path.getFileName().toString().toLowerCase(Locale.ROOT);
			if (relative.getNameCount() == 1 && Files.isDirectory(path, LinkOption.NOFOLLOW_LINKS)
					&& firstSegment(relative).equals("combined")) {
				return true;
			}
			if (relative.getNameCount() == 2 && Files.isRegularFile(path, LinkOption.NOFOLLOW_LINKS)
					&& name.equals(firstSegment(relative) + ".parquet")) {
				return true;
			}
			if (relative.getNameCount() == 2 && Files.isRegularFile(path, LinkOption.NOFOLLOW_LINKS)
					&& ((firstSegment(relative).equals("marker") && name.equals("marker.parquet"))
							|| (firstSegment(relative).equals("seed") && name.equals("seed.parquet")))) {
				return true;
			}
		}
		return false;
	}

	private static Set<String> allowedDiaMarkerRoots(String refinementMode) {
		String mode = refinementMode == null ? "" : refinementMode.toLowerCase(Locale.ROOT);
		if (mode.contains("deep")) {
			return Set.of("marker_deep");
		}
		if (mode.contains("comprehensive")) {
			return Set.of("marker_comprehensive");
		}
		if (mode.contains("focus") || mode.contains("core")) {
			return Set.of("marker_focus");
		}
		return Set.of("marker_focus", "marker_comprehensive", "marker_deep");
	}

	private static Set<String> allowedDiaSeedRoots(String refinementMode) {
		String mode = refinementMode == null ? "" : refinementMode.toLowerCase(Locale.ROOT);
		if (mode.contains("direct")) {
			return Set.of("seed_direct");
		}
		if (mode.contains("deep")) {
			return Set.of("seed_deep");
		}
		if (mode.contains("comprehensive")) {
			return Set.of("seed_comprehensive");
		}
		if (mode.contains("focus") || mode.contains("core")) {
			return Set.of("seed_focus");
		}
		return Set.of("seed_direct", "seed_focus", "seed_comprehensive", "seed_deep");
	}

	private static boolean containsSegment(Path relative, Set<String> names) {
		for (int i = 0; i < relative.getNameCount(); i++) {
			String segment = relative.getName(i).toString().toLowerCase(Locale.ROOT);
			if (names.contains(segment)) {
				return true;
			}
		}
		return false;
	}

	private static boolean isDiaWorkflow(WorkflowType workflowType) {
		return workflowType == WorkflowType.DIA_MARKER || workflowType == WorkflowType.DIA_SEED
				|| workflowType == WorkflowType.DIA_TARGET || workflowType == WorkflowType.DIA_REALTIME
				|| workflowType == WorkflowType.DIA_DEEP_COMPREHENSIVE;
	}

	private static WorkflowType workflowFromRootName(String rootName) {
		switch (rootName) {
		case "maxquant_search":
			return WorkflowType.DDA_IGC_MAXQUANT;
		case "pfind_igc_search":
			return WorkflowType.DDA_IGC_PFIND;
		case "fragpipe_igc_search":
			return WorkflowType.DDA_IGC_FRAGPIPE;
		case "sage_igc_search":
			return WorkflowType.DDA_IGC_SAGE;
		case "pfind_mag_search":
			return WorkflowType.DDA_MAG_PFIND;
		case "fragpipe_mag_search":
			return WorkflowType.DDA_MAG_FRAGPIPE;
		case "sage_mag_search":
			return WorkflowType.DDA_MAG_SAGE;
		case "marker_focus":
		case "marker_comprehensive":
		case "marker_deep":
			return WorkflowType.DIA_MARKER;
		case "seed_direct":
		case "seed_focus":
		case "seed_comprehensive":
		case "seed_deep":
			return WorkflowType.DIA_SEED;
		default:
			return null;
		}
	}

	private static String firstSegment(Path relative) {
		if (relative == null || relative.getNameCount() == 0) {
			return "";
		}
		return relative.getName(0).toString().toLowerCase(Locale.ROOT);
	}

	private static long safeSize(Path path) {
		try {
			return Files.size(path);
		} catch (IOException e) {
			return 0L;
		}
	}

	private static void addScore(Map<WorkflowType, Integer> scores, WorkflowType workflowType, int score) {
		scores.put(workflowType, scores.getOrDefault(workflowType, 0) + score);
	}

	private static JsonDetection detectFromJson(Path file) {
		try {
			String content = Files.readString(file, StandardCharsets.UTF_8);
			if (content.isBlank()) {
				return JsonDetection.empty();
			}
			JSONObject json = new JSONObject(content);
			String workflow = json.optString("workflow", null);
			WorkflowType workflowType = WorkflowType.fromAnalysisType(workflow);
			String refinementMode = firstPresent(json, "marker_refinement_mode", "seed_refinement_mode", "markerMode",
					"seedMode", "refinementMode");
			return new JsonDetection(workflowType, refinementMode);
		} catch (Exception e) {
			return JsonDetection.empty();
		}
	}

	private static String firstPresent(JSONObject json, String... keys) {
		for (String key : keys) {
			String value = json.optString(key, null);
			if (value != null && !value.isBlank()) {
				return value;
			}
		}
		return null;
	}

	private static Map<String, WorkflowType> createParameterFileWorkflowMap() {
		Map<String, WorkflowType> map = new LinkedHashMap<>();
		map.put("analysis_parameters.json", WorkflowType.DDA_IGC_MAXQUANT);
		map.put("analysis_parameters_pfind.json", WorkflowType.DDA_IGC_PFIND);
		map.put("analysis_parameters_fragpipe.json", WorkflowType.DDA_IGC_FRAGPIPE);
		map.put("analysis_parameters_sage.json", WorkflowType.DDA_IGC_SAGE);
		map.put("mag_pfind_analysis_parameters.json", WorkflowType.DDA_MAG_PFIND);
		map.put("mag_pfind_parameters.json", WorkflowType.DDA_MAG_PFIND);
		map.put("mag_fragpipe_analysis_parameters.json", WorkflowType.DDA_MAG_FRAGPIPE);
		map.put("mag_sage_analysis_parameters.json", WorkflowType.DDA_MAG_SAGE);
		map.put(uk.ac.quadram.metapilot.core.opensearch.OpenSearchLayout.PARAMETER_FILE,
				WorkflowType.DDA_MAG_OPENSEARCH);
		map.put("dia_marker_analysis_parameters.json", WorkflowType.DIA_MARKER);
		map.put("dia_seed_analysis_parameters.json", WorkflowType.DIA_SEED);
		map.put("dia_analysis_parameters.json", WorkflowType.DIA_SEED);
		map.put("dia_target_analysis_parameters.json", WorkflowType.DIA_TARGET);
		map.put("dia_rt_analysis_parameters.json", WorkflowType.DIA_REALTIME);
		return map;
	}

	private static final class JsonDetection {
		private final WorkflowType workflowType;
		private final String refinementMode;

		private JsonDetection(WorkflowType workflowType, String refinementMode) {
			this.workflowType = workflowType;
			this.refinementMode = refinementMode;
		}

		private static JsonDetection empty() {
			return new JsonDetection(null, null);
		}
	}

	private static void refuseUntilWorkflowProfilesVerified() throws IOException {
		if (!WORKFLOW_PROFILES_VERIFIED) {
			throw new IOException("Result folder cleaning is disabled because the current deletion rules are not "
					+ "workflow-verified. A safe cleaner must use explicit per-workflow keep/delete profiles derived "
					+ "from the runner code before preview or deletion is allowed.");
		}
	}

	private static final class KeepDecision {
		private final boolean keep;
		private final String reason;

		private KeepDecision(boolean keep, String reason) {
			this.keep = keep;
			this.reason = reason;
		}

		private static KeepDecision keep(String reason) {
			return new KeepDecision(true, reason);
		}

		private static KeepDecision delete() {
			return new KeepDecision(false, "Not retained by selected clean mode");
		}

		private static KeepDecision delete(String reason) {
			return new KeepDecision(false, reason);
		}
	}
}
