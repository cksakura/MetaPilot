package uk.ac.quadram.metapilot.task.dia;

import java.io.File;
import java.util.ArrayList;

import javafx.application.Platform;
import javafx.scene.control.ProgressBar;
import uk.ac.quadram.metapilot.resources.MetaResourceIO;
import uk.ac.quadram.metapilot.resources.MetaResourceParameters;
import uk.ac.quadram.metapilot.task.dia.par.MetaParameterDiaMarker;
import uk.ac.quadram.metapilot.task.mag.MagDbConfigIO;
import uk.ac.quadram.metapilot.task.mag.MagDbItem;
import uk.ac.quadram.metapilot.task.mag.MagDbItem.MagVersion;
import uk.ac.quadram.metapilot.task.par.MetaBaseParameter;
import uk.ac.quadram.metapilot.task.par.MetaData;
import uk.ac.quadram.metapilot.workflow.WorkflowType;

/**
 * Minimal benchmark harness for testing deep marker strategies on a single raw
 * file.
 *
 * Intended usage: - Run from IDE with main() or call runOnce(...) - Provide raw
 * file, output folder, MAG DB config (MAGDB_1_0.json), catalogue ID/version -
 * Optionally set genomesPerRound (e.g. 5 vs 10)
 *
 * This runner keeps parameters simple: defaults everywhere except what you
 * pass.
 */
public final class DeepMarkerSingleRawBenchmark {

	private DeepMarkerSingleRawBenchmark() {
	}

	private static volatile boolean fxStarted = false;

	private static synchronized void ensureJavaFxStarted() {
		if (fxStarted)
			return;
		try {
			Platform.startup(() -> {
				// no-op: just initialize toolkit
			});
		} catch (IllegalStateException alreadyStarted) {
			// Toolkit already initialized in this JVM
		}
		fxStarted = true;
	}

	public static void main(String[] args) throws Exception {
		if (args.length < 5) {
			System.err.println("Usage:");
			System.err.println("  <rawFile> <outDir> <magdbJson> <catalogueId> <catalogueVersion> [genomesPerRound]");
			System.err.println();
			System.err.println("Example:");
			System.err.println("  D:\\data\\run1.raw  D:\\out  C:\\...\\MAGDB_1_0.json  MGnify  v1  10");
			System.exit(2);
		}
		String rawFile = args[0];
		String outDir = args[1];
		String magdbJson = args[2];
		String catalogueId = args[3];
		String catalogueVersion = args[4];
		int genomesPerRound = args.length >= 6 ? Integer.parseInt(args[5]) : 5;

		runOnce(rawFile, outDir, magdbJson, catalogueId, catalogueVersion, genomesPerRound);
	}

	public static boolean runOnce(String rawFilePath, String outDirPath, String magDbConfigJsonPath,
			String magDbCatalogueId, String magDbCatalogueVersion, int genomesPerRound) throws Exception {
		ensureJavaFxStarted();

		File raw = new File(rawFilePath);
		File outDir = new File(outDirPath);
		File magDbJson = new File(magDbConfigJsonPath);
		if (!raw.exists())
			throw new IllegalArgumentException("Raw file not found: " + raw.getAbsolutePath());
		if (!magDbJson.exists())
			throw new IllegalArgumentException("MAG DB config not found: " + magDbJson.getAbsolutePath());
		outDir.mkdirs();

		MagDbItem magDb = resolveMagDb(magDbJson, magDbCatalogueId, magDbCatalogueVersion);
		if (magDb == null)
			throw new IllegalStateException("Failed to resolve MAG DB item from config");

		MetaResourceParameters mrp = new MetaResourceIO().loadParameters();

		MetaData md = new MetaData(new String[] { raw.getAbsolutePath() }, new String[] { raw.getName() });
		MetaBaseParameter base = new MetaBaseParameter(WorkflowType.DIA_MARKER.getAnalysisType(), md,
				outDir.getAbsolutePath(), magDbCatalogueId, false, "", "DIA", 1,
				Math.max(1, Runtime.getRuntime().availableProcessors() / 2), true);

		MetaParameterDiaMarker p = new MetaParameterDiaMarker(base, magDb, true, false, new String[0], false,
				new String[0]);
		p.setMarkerMode("DEEP_PER_FILE");
		p.setGenomesPerRound(genomesPerRound <= 0 ? 5 : genomesPerRound);
		long t0 = System.currentTimeMillis();
		DiaMarkerDeepPerFileTask task = new DiaMarkerDeepPerFileTask(p, mrp, new ProgressBar(), new ProgressBar());
		boolean ok = Boolean.TRUE.equals(task.call());
		long t1 = System.currentTimeMillis();

		System.out.println("DEEP_PER_FILE single-raw benchmark finished.");
		System.out.println("  ok=" + ok);
		System.out.println("  genomesPerRound=" + p.getGenomesPerRound());
		System.out.println("  elapsed_ms=" + (t1 - t0));
		System.out.println("  outDir=" + outDir.getAbsolutePath());
		return ok;
	}

	private static MagDbItem resolveMagDb(File magDbJson, String catalogueId, String catalogueVersion) {
		MagDbItem[][] items = MagDbConfigIO.parse(magDbJson);
		if (items == null)
			return null;

		ArrayList<MagDbItem> flat = new ArrayList<>();
		for (MagDbItem[] group : items) {
			if (group == null)
				continue;
			for (MagDbItem it : group) {
				if (it != null)
					flat.add(it);
			}
		}
		for (MagDbItem it : flat) {
			if (it.getCatalogueID() == null)
				continue;
			if (!it.getCatalogueID().equals(catalogueId))
				continue;
			if (it.getVersionList() == null)
				continue;
			for (MagVersion v : it.getVersionList()) {
				if (v == null)
					continue;
				if (v.getCatalogueVersion() != null && v.getCatalogueVersion().equals(catalogueVersion)) {
					it.setUsedVersion(v);
					return it;
				}
			}
		}
		return null;
	}
}
