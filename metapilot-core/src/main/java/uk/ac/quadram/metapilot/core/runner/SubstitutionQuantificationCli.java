package uk.ac.quadram.metapilot.core.runner;

import java.io.File;
import java.nio.file.Paths;

import uk.ac.quadram.metapilot.core.resources.MetaPilotResourceConfig;
import uk.ac.quadram.metapilot.core.resources.ToolKey;
import uk.ac.quadram.metapilot.core.workflow.ProgressCallback;
import uk.ac.quadram.metapilot.task.dda.mag.par.DDAMagOpenSearchParameterIO;
import uk.ac.quadram.metapilot.task.dda.mag.par.DDAMagOpenSearchParameters;

/**
 * Runs only the side quantification in {@code sequence_variants/} on a finished open-search result
 * ({@link SubstitutionQuantification}): no search, identification or main quantification runs, and
 * no main table is written. For checking the side run on existing results.
 *
 * <pre>
 * java -cp metapilot-cli-shaded.jar uk.ac.quadram.metapilot.core.runner.SubstitutionQuantificationCli \
 *     &lt;opensearch_MAG_search folder&gt; --params workflow.json --resources meta_resources.params [--threads N]
 * </pre>
 */
public final class SubstitutionQuantificationCli {

	private SubstitutionQuantificationCli() {
	}

	public static void main(String[] args) throws Exception {
		if (args.length < 5) {
			System.err.println("Usage: SubstitutionQuantificationCli <opensearch_MAG_search folder> --params workflow.json "
					+ "--resources meta_resources.params [--threads N]");
			System.exit(2);
		}
		File root = new File(args[0]);
		File paramsFile = null;
		File resources = null;
		int threads = Math.max(1, Runtime.getRuntime().availableProcessors());
		for (int i = 1; i < args.length; i++) {
			if ("--params".equals(args[i]) && i + 1 < args.length) {
				paramsFile = new File(args[++i]);
			} else if ("--resources".equals(args[i]) && i + 1 < args.length) {
				resources = new File(args[++i]);
			} else if ("--threads".equals(args[i]) && i + 1 < args.length) {
				threads = Integer.parseInt(args[++i]);
			} else {
				System.err.println("Unknown argument: " + args[i]);
				System.exit(2);
			}
		}
		if (paramsFile == null || resources == null || !root.isDirectory()) {
			System.err.println("--params and --resources are required, and the result folder must exist: " + root);
			System.exit(2);
		}
		DDAMagOpenSearchParameters params = new DDAMagOpenSearchParameterIO().loadParameters(paramsFile);
		File flashLfq = MetaPilotResourceConfig.load(Paths.get(resources.getPath()))
				.resolveToolPath(ToolKey.FLASHLFQ, null).map(java.nio.file.Path::toFile).orElse(null);
		System.out.println("FlashLFQ: " + flashLfq);
		ProgressCallback console = new ProgressCallback() {
			@Override
			public void onProgress(double workDone, double max) {
			}

			@Override
			public void onMessage(String message) {
				System.out.println(message);
			}

			@Override
			public void onComplete(Object result) {
			}

			@Override
			public void onError(Exception error) {
				error.printStackTrace(System.err);
			}
		};
		File spectraDir = MagOpenSearchHeadlessRunner.readableSpectraDir(params.getExperimentData(), root);
		System.out.println("Spectra: " + spectraDir);
		File out = MagOpenSearchHeadlessRunner.quantifySubstitutionPsms(params, flashLfq, threads, root,
				params.getSearchEngines(), spectraDir, console);
		if (out == null) {
			System.err.println("No side quantification table was written (see the messages above)");
			System.exit(1);
		}
		System.out.println("Side quantification: " + out.getAbsolutePath());
	}
}
