package uk.ac.quadram.metapilot.core.runner;

import java.io.File;
import java.util.ArrayList;
import java.util.List;

import uk.ac.quadram.metapilot.core.opensearch.ModificationValidator;
import uk.ac.quadram.metapilot.core.workflow.ProgressCallback;
import uk.ac.quadram.metapilot.task.dda.mag.par.DDAMagOpenSearchParameterIO;
import uk.ac.quadram.metapilot.task.dda.mag.par.DDAMagOpenSearchParameters;

/** Command-line entry point for rerunning pooled validation from existing pepXML. */
public final class OpenSearchValidationCli {

	private OpenSearchValidationCli() {
	}

	public static void main(String[] args) throws Exception {
		if (args.length < 2) {
			throw new IllegalArgumentException(
					"Usage: OpenSearchValidationCli <output-root> [--params workflow.json] <file1.pepXML> ...");
		}
		File outputRoot = new File(args[0]);
		DDAMagOpenSearchParameters params = new DDAMagOpenSearchParameters();
		List<File> pepXmlFiles = new ArrayList<>();
		List<File> sageResults = new ArrayList<>();
		for (int i = 1; i < args.length; i++) {
			if ("--params".equals(args[i])) {
				if (++i >= args.length) {
					throw new IllegalArgumentException("--params requires a workflow JSON path");
				}
				params = new DDAMagOpenSearchParameterIO().loadParameters(new File(args[i]));
			} else if ("--sage".equals(args[i]) || "--engines".equals(args[i])) {
				throw new IllegalArgumentException("MSFragger + Sage has been removed. Closed and open "
						+ "MSFragger searches replace a second engine; do not pass --sage or --engines.");
			} else {
				pepXmlFiles.add(new File(args[i]));
			}
		}

		ProgressCallback console = new ProgressCallback() {
			@Override
			public void onProgress(double workDone, double max) {
				// The validation internals report meaningful textual milestones.
			}

			@Override
			public void onMessage(String message) {
				System.out.println(message);
			}

			@Override
			public void onComplete(Object result) {
				System.out.println("Complete: " + result);
			}

			@Override
			public void onError(Exception error) {
				error.printStackTrace(System.err);
			}
		};

		params.setSearchEngines(DDAMagOpenSearchParameters.ENGINES_MSFRAGGER);
		System.out.println("Engine: MSFragger; " + pepXmlFiles.size() + " pepXML file(s)");

		ModificationValidator.Result result = MagOpenSearchHeadlessRunner.validatePepXmls(params, pepXmlFiles,
				sageResults, outputRoot, console);
		System.out.println("Validation complete: " + result.validated.size() + " retained PSMs, "
				+ result.acceptedPeakCount() + " accepted fitted peaks; outputs: "
				+ new File(outputRoot, "combined").getAbsolutePath());
	}
}
