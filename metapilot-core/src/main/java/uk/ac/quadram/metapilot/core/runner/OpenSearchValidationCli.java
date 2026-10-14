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
					"Usage: OpenSearchValidationCli <output-root> [--params workflow.json] "
							+ "[--engines \"MSFragger + Sage\"] [--sage results.sage.tsv ...] <file1.pepXML> ...");
		}
		File outputRoot = new File(args[0]);
		DDAMagOpenSearchParameters params = new DDAMagOpenSearchParameters();
		List<File> pepXmlFiles = new ArrayList<>();
		// Sage results are accepted here so a MSFragger+Sage run can be revalidated
		// from existing searches. Without this the CLI could only ever reproduce a
		// MSFragger-only result, which is not the configuration that needs the
		// evaluation loop most.
		List<File> sageResults = new ArrayList<>();
		String engines = null;
		for (int i = 1; i < args.length; i++) {
			if ("--params".equals(args[i])) {
				if (++i >= args.length) {
					throw new IllegalArgumentException("--params requires a workflow JSON path");
				}
				params = new DDAMagOpenSearchParameterIO().loadParameters(new File(args[i]));
			} else if ("--sage".equals(args[i])) {
				if (++i >= args.length) {
					throw new IllegalArgumentException("--sage requires a results.sage.tsv path");
				}
				sageResults.add(new File(args[i]));
			} else if ("--engines".equals(args[i])) {
				if (++i >= args.length) {
					throw new IllegalArgumentException("--engines requires an engine-mode string");
				}
				engines = args[i];
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

		// Engine mode decides which combined_* folder the outputs land in, so it
		// must match the run being reproduced or the results go somewhere else.
		if (engines != null) {
			params.setSearchEngines(engines);
		} else if (!sageResults.isEmpty()) {
			params.setSearchEngines("MSFragger + Sage");
		}
		System.out.println("Engine mode: " + params.getSearchEngines() + "; " + pepXmlFiles.size()
				+ " pepXML file(s), " + sageResults.size() + " Sage result file(s)");

		ModificationValidator.Result result = MagOpenSearchHeadlessRunner.validatePepXmls(params, pepXmlFiles,
				sageResults, outputRoot, console);
		System.out.println("Validation complete: " + result.validated.size() + " retained PSMs, "
				+ result.acceptedPeakCount() + " accepted fitted peaks; outputs: "
				+ new File(outputRoot, "combined").getAbsolutePath());
	}
}
