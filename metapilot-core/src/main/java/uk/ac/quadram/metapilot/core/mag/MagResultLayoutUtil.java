package uk.ac.quadram.metapilot.core.mag;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.StandardCopyOption;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashSet;
import java.util.Set;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

/**
 * Shared helpers for the "nested search_results" DDA MAG output layout.
 *
 * <p>
 * Historically every DDA MAG engine (pFind, FragPipe, Sage) wrote one folder
 * per raw/sample file directly under the workflow result root, as a sibling
 * of the fixed {@code combined/}, {@code taxonomy_analysis/},
 * {@code functional_annotation/}, and {@code report/} folders. With many raw
 * files this makes the fixed, always-useful folders hard to find among the
 * per-sample ones - the DIA workflows solved the same problem by nesting all
 * per-file folders one level deeper (see e.g.
 * {@code DiaMarkerDeepPerFileHeadlessRunner}, which nests per-file folders
 * under {@code marker/}).
 *
 * <p>
 * This class introduces the equivalent nested layout for DDA MAG, opt-in per
 * engine via {@code MetaPilotMagTask.useNestedSearchResultsLayout()}:
 *
 * <pre>
 * &lt;result&gt;/&lt;engine&gt;_MAG_search/
 *   search_results/&lt;rawFileName&gt;/...   &lt;- per-raw-file search output, nested
 *   combined/                            &lt;- FlashLFQ / quant combine output
 *   taxonomy_analysis/
 *   functional_annotation/
 *   excluded_evidence/
 *   report/
 * </pre>
 *
 * <p>
 * All DDA MAG engines have opted in: pFind ({@code DDAMagPfindTask}),
 * Sage ({@code DDAMagSageTask}), FragPipe ({@code MagFragPipeHeadlessRunner}),
 * and open search ({@code MagOpenSearchHeadlessRunner} /
 * {@code DDAMagOpenSearchTask}).
 *
 * <p>
 * <b>Wrapper folder name:</b> deliberately {@code search_results}, not
 * {@code marker} - each per-sample folder already contains its own
 * {@code marker/} subfolder (marker-gene search results), so reusing
 * {@code marker} as the top-level wrapper would produce a confusing
 * {@code search_results/<sample>/marker/} vs. a hypothetical
 * {@code marker/<sample>/marker/} naming collision.
 */
public final class MagResultLayoutUtil {

	private static final Logger logger = LogManager.getLogger(MagResultLayoutUtil.class);

	public static final String SEARCH_RESULTS_DIR_NAME = "search_results";

	/**
	 * pFind writes {@code ms2Count.txt} and {@code rt.txt} in the result's
	 * {@code spectra/} folder. Open search and FragPipe do the same. This name
	 * is only the earlier location under the engine root; a leftover copy is
	 * still read, and it stays in {@link #FIXED_TOP_LEVEL_DIR_NAMES} so the
	 * nested-layout migration does not hide it under {@code search_results/}.
	 */
	public static final String SPECTRA_SUMMARY_DIR_NAME = "spectra_summary";

	/** {@code <result>/spectra}, where pFind keeps the count files. */
	public static final String SPECTRA_DIR_NAME = "spectra";

	/**
	 * Peptides dropped after the genome set is frozen. Written beside
	 * {@code taxonomy_analysis/} by every MAG workflow, and kept there so the
	 * nested-layout migration does not treat it as a sample folder.
	 */
	public static final String EXCLUDED_EVIDENCE_DIR_NAME = "excluded_evidence";

	private static final Set<String> FIXED_TOP_LEVEL_DIR_NAMES = Collections
			.unmodifiableSet(new HashSet<>(Arrays.asList("combined", "taxonomy_analysis", "functional_annotation",
					"report", "host", SEARCH_RESULTS_DIR_NAME, SPECTRA_SUMMARY_DIR_NAME,
					EXCLUDED_EVIDENCE_DIR_NAME)));

	private MagResultLayoutUtil() {
	}

	/**
	 * Returns the directory that per-raw-file sample folders should live under
	 * for engines using the nested layout, migrating any legacy flat-layout
	 * sample folders into it first - regardless of whether that prior run
	 * finished quantification (see
	 * {@link #migrateLegacyFlatLayoutIfNeeded(File, File)}) - then creating it
	 * if it still does not exist.
	 *
	 * @param resultRoot e.g. {@code <result>/pfind_MAG_search}
	 */
	public static File resolveSearchResultsRoot(File resultRoot) {
		return resolveSearchResultsRoot(resultRoot, Collections.emptySet());
	}

	/**
	 * As {@link #resolveSearchResultsRoot(File)}, leaving the workflow's own
	 * top-level folders {@code workflowDirNames} where they are rather than
	 * taking them for legacy sample folders.
	 */
	public static File resolveSearchResultsRoot(File resultRoot, Set<String> workflowDirNames) {
		File searchResultsDir = new File(resultRoot, SEARCH_RESULTS_DIR_NAME);
		migrateLegacyFlatLayoutIfNeeded(resultRoot, searchResultsDir, workflowDirNames);
		if (!searchResultsDir.exists()) {
			searchResultsDir.mkdirs();
		}
		return searchResultsDir;
	}

	/**
	 * Moves legacy per-raw-file folders that sit directly under
	 * {@code resultRoot} (from before an engine opted into the nested layout)
	 * into {@code searchResultsDir}.
	 *
	 * <p>
	 * Unconditional on completion of the run - deliberately NOT gated on
	 * {@code combined/QuantifiedPeptides.tsv}/{@code psms.tsv} existing.
	 * {@code DDAMagPfindTask} has several per-sample resume checks that run
	 * <em>before</em> quantification (marker search, sample-specific FASTA
	 * extraction, sample-DB search - see {@code separateSearchMarker()},
	 * {@code fastaExtractWithGenomeSelectorMode()}, {@code separateSearch()}),
	 * so a legacy flat-layout run interrupted between "search done" and
	 * "quant done" must still be migrated, or those already-completed
	 * per-sample steps would be silently re-run against an empty
	 * {@code search_results/}. FragPipe MAG still resolves this once at the
	 * start of {@code runInternal()}. Open search calls it throughout the run
	 * via {@code OpenSearchLayout.searchResultsRoot}, so names in
	 * {@link #FIXED_TOP_LEVEL_DIR_NAMES} (including {@code spectra_summary/})
	 * must stay here or they are moved as leftover sample folders. Running two
	 * task instances against the same result folder at once is unsupported
	 * regardless of this method.
	 */
	public static void migrateLegacyFlatLayoutIfNeeded(File resultRoot, File searchResultsDir) {
		migrateLegacyFlatLayoutIfNeeded(resultRoot, searchResultsDir, Collections.emptySet());
	}

	private static void migrateLegacyFlatLayoutIfNeeded(File resultRoot, File searchResultsDir,
			Set<String> workflowDirNames) {
		if (resultRoot == null || !resultRoot.isDirectory()) {
			return;
		}

		File[] children = resultRoot.listFiles(File::isDirectory);
		if (children == null) {
			return;
		}

		for (File child : children) {
			String name = child.getName();
			if (FIXED_TOP_LEVEL_DIR_NAMES.contains(name) || workflowDirNames.contains(name)) {
				continue;
			}

			// A directory directly under resultRoot that isn't one of the fixed names -
			// this is a legacy per-raw-file sample folder. Move it under search_results/.
			if (!searchResultsDir.exists()) {
				searchResultsDir.mkdirs();
			}
			File destination = new File(searchResultsDir, name);
			if (destination.exists()) {
				logger.warn("Migration skipped for {} - destination {} already exists", child.getAbsolutePath(),
						destination.getAbsolutePath());
				continue;
			}
			try {
				Files.move(child.toPath(), destination.toPath(), StandardCopyOption.ATOMIC_MOVE);
				logger.info("Migrated legacy DDA MAG result folder {} -> {}", child.getAbsolutePath(),
						destination.getAbsolutePath());
			} catch (IOException atomicMoveFailed) {
				try {
					Files.move(child.toPath(), destination.toPath());
					logger.info("Migrated legacy DDA MAG result folder {} -> {} (non-atomic fallback)",
							child.getAbsolutePath(), destination.getAbsolutePath());
				} catch (IOException e2) {
					logger.error("Failed to migrate legacy DDA MAG result folder {} -> {}: {}",
							child.getAbsolutePath(), destination.getAbsolutePath(), e2.getMessage(), e2);
				}
			}
		}
	}

	/**
	 * Names that are never treated as a per-sample folder when scanning a DDA
	 * MAG result root (flat legacy layout) or a {@code search_results/} folder.
	 */
	public static Set<String> fixedTopLevelDirNames() {
		return FIXED_TOP_LEVEL_DIR_NAMES;
	}

	/**
	 * {@code <result>/spectra}, the parent of an engine folder such as
	 * {@code opensearch_MAG_search}. This is {@code MetaBaseParameter.getSpectraFile()}.
	 */
	public static File spectraDir(File engineResultRoot) {
		if (engineResultRoot == null) {
			return null;
		}
		File result = engineResultRoot.getParentFile();
		return new File(result == null ? engineResultRoot : result, SPECTRA_DIR_NAME);
	}

	/** Earlier location, {@code <engineResultRoot>/spectra_summary}. */
	public static File spectraSummaryDir(File resultRoot) {
		return new File(resultRoot, SPECTRA_SUMMARY_DIR_NAME);
	}

	public static boolean hasMs2Count(File dir) {
		if (dir == null) {
			return false;
		}
		File file = new File(dir, "ms2Count.txt");
		return file.isFile() && file.length() > 0;
	}

	/**
	 * {@code ms2Count.txt} for the identification rate. Checks the folder itself,
	 * {@code spectra/} inside it, {@code spectra/} beside an engine folder
	 * (pFind's place, and where new open-search and FragPipe runs write), then
	 * a leftover {@code spectra_summary/} at the engine root or under
	 * {@code search_results/}.
	 */
	public static File findMs2Count(File resultRoot) {
		return locateMs2Count(resultRoot, null);
	}

	/**
	 * Same search as {@link #findMs2Count(File)}, starting from the spectra
	 * folder the caller already has and then the engine result folder.
	 */
	public static File locateMs2Count(File spectraDir, File engineResultDir) {
		File found = countNear(spectraDir);
		if (found != null) {
			return found;
		}
		return countNear(engineResultDir);
	}

	private static File countNear(File dir) {
		if (dir == null) {
			return null;
		}
		File direct = countFile(dir);
		if (direct != null) {
			return direct;
		}
		File inside = countFile(new File(dir, SPECTRA_DIR_NAME));
		if (inside != null) {
			return inside;
		}
		File besideEngine = countFile(spectraDir(dir));
		if (besideEngine != null) {
			return besideEngine;
		}
		File legacy = countFile(spectraSummaryDir(dir));
		if (legacy != null) {
			return legacy;
		}
		if (SPECTRA_SUMMARY_DIR_NAME.equals(dir.getName()) && dir.getParentFile() != null) {
			File parent = dir.getParentFile();
			File engine = SEARCH_RESULTS_DIR_NAME.equals(parent.getName()) ? parent.getParentFile() : parent;
			File fromSummaryParent = countNearEngineOnly(engine);
			if (fromSummaryParent != null) {
				return fromSummaryParent;
			}
		}
		return countFile(new File(new File(dir, SEARCH_RESULTS_DIR_NAME), SPECTRA_SUMMARY_DIR_NAME));
	}

	/** Engine-root locations only, so a lookup that started in {@code spectra_summary/} cannot recurse. */
	private static File countNearEngineOnly(File engine) {
		if (engine == null) {
			return null;
		}
		File beside = countFile(spectraDir(engine));
		if (beside != null) {
			return beside;
		}
		File legacy = countFile(spectraSummaryDir(engine));
		if (legacy != null) {
			return legacy;
		}
		return countFile(new File(new File(engine, SEARCH_RESULTS_DIR_NAME), SPECTRA_SUMMARY_DIR_NAME));
	}

	private static File countFile(File dir) {
		if (!hasMs2Count(dir)) {
			return null;
		}
		return new File(dir, "ms2Count.txt");
	}

	/**
	 * The folder that should hold {@code ms2Count.txt} and {@code rt.txt}:
	 * {@code <result>/spectra}, matching pFind. A copy already in
	 * {@code spectra_summary/} (including one nested under {@code search_results/})
	 * is moved there when {@code spectra/} does not have the count yet.
	 */
	public static File spectraCountDir(File engineResultRoot) {
		File spectra = spectraDir(engineResultRoot);
		if (spectra == null || hasMs2Count(spectra)) {
			return spectra;
		}
		File legacy = legacyCountDir(engineResultRoot);
		if (legacy == null) {
			return spectra;
		}
		if (!spectra.isDirectory() && !spectra.mkdirs()) {
			logger.warn("Could not create {}; using {}", spectra.getAbsolutePath(), legacy.getAbsolutePath());
			return legacy;
		}
		boolean moved = moveCountFile(legacy, spectra, "ms2Count.txt") && moveCountFile(legacy, spectra, "rt.txt");
		if (moved || hasMs2Count(spectra)) {
			logger.info("Moved spectrum counts {} -> {}", legacy.getAbsolutePath(), spectra.getAbsolutePath());
			return spectra;
		}
		return legacy;
	}

	private static File legacyCountDir(File engineResultRoot) {
		if (engineResultRoot == null) {
			return null;
		}
		File atRoot = spectraSummaryDir(engineResultRoot);
		if (hasMs2Count(atRoot)) {
			return atRoot;
		}
		File nested = new File(new File(engineResultRoot, SEARCH_RESULTS_DIR_NAME), SPECTRA_SUMMARY_DIR_NAME);
		return hasMs2Count(nested) ? nested : null;
	}

	/** True when the file was absent at the source, or now sits in {@code destDir}. */
	private static boolean moveCountFile(File fromDir, File destDir, String name) {
		File from = new File(fromDir, name);
		if (!from.isFile()) {
			return true;
		}
		File dest = new File(destDir, name);
		if (dest.isFile()) {
			return true;
		}
		try {
			Files.move(from.toPath(), dest.toPath(), StandardCopyOption.ATOMIC_MOVE);
			return true;
		} catch (IOException atomicMoveFailed) {
			try {
				Files.move(from.toPath(), dest.toPath());
				return true;
			} catch (IOException e2) {
				logger.warn("Could not move {} -> {}: {}", from.getAbsolutePath(), dest.getAbsolutePath(),
						e2.getMessage());
				return false;
			}
		}
	}

}
