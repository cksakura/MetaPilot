package uk.ac.quadram.metapilot.core.unified;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

import uk.ac.quadram.metapilot.core.opensearch.OpenSearchPepxmlReader;
import uk.ac.quadram.metapilot.core.opensearch.OpenSearchPsm;
import uk.ac.quadram.metapilot.core.prodb.TargetDecoyFasta;
import uk.ac.quadram.metapilot.task.dda.mag.par.DDAMagOpenSearchParameterIO;
import uk.ac.quadram.metapilot.task.dda.mag.par.DDAMagOpenSearchParameters;

/**
 * Runs the unified identification on existing searches, writing only into {@code --out}.
 *
 * <pre>
 * UnifiedIdentificationCli --params workflow.json --open open1.pepXML,... --closed closed1.pepXML,...
 *     --mzml dir --out dir [--diann DiaNN.exe] [--closed-topn 3] [--report-unexplained]
 *     [--max-training-queries N] [--stagewise] [--psm-recurrence]
 *     diagnostics: [--closed-only] [--engine-features] [--no-reset]
 * </pre>
 *
 * {@code --open}/{@code --closed} also accept a directory (every *.pepXML inside).
 */
public final class UnifiedIdentificationCli {

	private UnifiedIdentificationCli() {
	}

	public static void main(String[] args) throws Exception {
		Map<String, String> a = new LinkedHashMap<>();
		for (int i = 0; i < args.length; i++) {
			if (args[i].equals("--report-unexplained") || args[i].equals("--closed-only")
					|| args[i].equals("--engine-features") || args[i].equals("--no-reset")
					|| args[i].equals("--stagewise") || args[i].equals("--psm-recurrence")) {
				a.put(args[i].substring(2), "true");
			} else if (args[i].startsWith("--") && i + 1 < args.length) {
				a.put(args[i].substring(2), args[++i]);
			}
		}
		UnifiedIdentification.Config cfg = new UnifiedIdentification.Config();
		if (a.containsKey("params")) {
			DDAMagOpenSearchParameters params = new DDAMagOpenSearchParameterIO().loadParameters(new File(a.get("params")));
			cfg.fitter = params.toFitterConfig();
			cfg.folds = params.getClassifierFolds();
		}
		cfg.reportUnexplained = a.containsKey("report-unexplained");
		cfg.stagewise = a.containsKey("stagewise");
		cfg.psmRecurrence = a.containsKey("psm-recurrence");
		// Diagnostics for checking the machinery against the closed classifier; not the method.
		cfg.closedOnly = a.containsKey("closed-only");
		cfg.engineFeatures = a.containsKey("engine-features");
		cfg.reset = !a.containsKey("no-reset");
		if (a.containsKey("drop-features")) {
			cfg.dropFeatures.addAll(Arrays.asList(a.get("drop-features").split(",")));
		}
		if (a.containsKey("max-training-queries")) {
			cfg.maxTrainingQueries = Integer.parseInt(a.get("max-training-queries"));
		}
		File out = new File(a.getOrDefault("out", "unified_out"));
		cfg.workDir = new File(out, "work");
		if (a.containsKey("diann")) {
			cfg.diann = new File(a.get("diann"));
		}
		if (a.containsKey("threads")) {
			cfg.threads = Integer.parseInt(a.get("threads"));
		}

		OpenSearchPepxmlReader openReader = new OpenSearchPepxmlReader(TargetDecoyFasta.DECOY_TAG);
		List<OpenSearchPsm> open = openReader.read(files(a.get("open")), System.out::println);
		OpenSearchPepxmlReader closedReader = new OpenSearchPepxmlReader(TargetDecoyFasta.DECOY_TAG);
		closedReader.setClosedSearch(true);
		closedReader.setMaxHitRank(Integer.parseInt(a.getOrDefault("closed-topn", "3")));
		List<OpenSearchPsm> closed = closedReader.read(files(a.get("closed")), System.out::println);

		Map<String, File> mzml = new LinkedHashMap<>();
		File[] spectra = new File(a.get("mzml")).listFiles((d, n) -> n.toLowerCase().endsWith(".mzml"));
		if (spectra != null) {
			for (File f : spectra) {
				mzml.put(f.getName().substring(0, f.getName().length() - 5), f);
			}
		}
		long t0 = System.currentTimeMillis();
		UnifiedIdentification.run(open, closed, mzml, cfg, out, System.out::println);
		Runtime rt = Runtime.getRuntime();
		System.out.printf("Unified DONE in %.1f min; heap used %.1f GB of %.1f GB; outputs in %s%n",
				(System.currentTimeMillis() - t0) / 60000.0, (rt.totalMemory() - rt.freeMemory()) / 1e9,
				rt.maxMemory() / 1e9, out.getAbsolutePath());
	}

	private static List<File> files(String spec) {
		List<File> out = new ArrayList<>();
		for (String part : spec.split(",")) {
			File f = new File(part);
			if (f.isDirectory()) {
				File[] inside = f.listFiles((d, n) -> n.endsWith(".pepXML"));
				if (inside != null) {
					Arrays.sort(inside);
					out.addAll(Arrays.asList(inside));
				}
			} else {
				out.add(f);
			}
		}
		return out;
	}
}
