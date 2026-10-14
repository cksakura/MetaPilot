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
 *   report/
 * </pre>
 *
 * <p>
 * All three DDA MAG engines have opted in: pFind ({@code DDAMagPfindTask}),
 * Sage ({@code DDAMagSageTask}), and FragPipe ({@code MagFragPipeHeadlessRunner},
 * which is also the single execution path used by {@code DDAMagFragPipeTask}
 * and every CLI/HPC FragPipe caller).
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

	private static final Set<String> FIXED_TOP_LEVEL_DIR_NAMES = Collections
			.unmodifiableSet(new HashSet<>(Arrays.asList("combined", "taxonomy_analysis", "functional_annotation",
					"report", "host", SEARCH_RESULTS_DIR_NAME)));

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
	 * {@code search_results/}. Safety instead comes from when this runs: it is
	 * only ever called once, synchronously, from
	 * {@code MetaPilotMagTask.createSeparateFolders()} at task construction -
	 * before any search step starts - so there is no concurrent writer for
	 * this task's own run to race. (Running two task instances against the
	 * same result folder at once is unsupported regardless of this method.)
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
}
