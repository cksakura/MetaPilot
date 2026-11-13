package uk.ac.quadram.metapilot.core.runner;

import uk.ac.quadram.metapilot.core.resources.MetaPilotResourceConfig;
import uk.ac.quadram.metapilot.core.resources.ToolKey;
import uk.ac.quadram.metapilot.resources.MetaResourceParameters;
import uk.ac.quadram.metapilot.resources.MetaResourceParameters.ToolType;
import uk.ac.quadram.metapilot.task.dda.mag.par.DDAMagOpenSearchParameters;
import uk.ac.quadram.metapilot.task.mag.MagDbItem;

import java.io.File;
import java.io.IOException;
import java.nio.file.Path;

/**
 * Resolves the external tools the MAG open-search workflow needs from the same
 * Resources entries the Windows GUI uses ({@link MetaResourceParameters}), then
 * writes the resolved paths onto the parameter object.
 * <p>
 * Shared by the GUI task, the {@code mag-opensearch} CLI command and the HPC
 * study runner so all three fail for the same reason when a tool is missing.
 */
public final class MagOpenSearchToolSetup {

	private MagOpenSearchToolSetup() {
	}

	public static final class Tools {
		public final File msfraggerJar;
		public final File markerFasta;
		public final File flashLfqExe;
		public final File funcDbFile;

		Tools(File msfraggerJar, File markerFasta, File flashLfqExe, File funcDbFile) {
			this.msfraggerJar = msfraggerJar;
			this.markerFasta = markerFasta;
			this.flashLfqExe = flashLfqExe;
			this.funcDbFile = funcDbFile;
		}
	}

	/**
	 * Builds a {@link MetaResourceParameters} from the CLI/HPC site resources
	 * file so the same FragPipe-bundled jar lookup the GUI uses also works
	 * headless.
	 */
	public static MetaResourceParameters fromResourceConfig(MetaPilotResourceConfig rc) {
		MetaResourceParameters mrp = new MetaResourceParameters();
		if (rc == null) {
			return mrp;
		}
		copy(rc, ToolKey.FRAGPIPE, ToolType.FRAGPIPE, mrp);
		copy(rc, ToolKey.FLASHLFQ, ToolType.FLASHLFQ, mrp);
		copy(rc, ToolKey.FUNCDB, ToolType.FUNCDB, mrp);
		// Optional keys may appear in meta_resources.params even when not in ToolKey.
		rc.getRaw("SAGE").filter(s -> s != null && !s.isBlank())
				.ifPresent(s -> mrp.setToolPath(ToolType.SAGE, s.trim()));
		rc.getRaw("MSCONVERT").filter(s -> s != null && !s.isBlank())
				.ifPresent(s -> mrp.setToolPath(ToolType.MSCONVERT, s.trim()));
		return mrp;
	}

	private static void copy(MetaPilotResourceConfig rc, ToolKey key, ToolType type, MetaResourceParameters mrp) {
		rc.resolveToolPath(key, null).map(Path::toString).ifPresent(path -> mrp.setToolPath(type, path));
	}

	/**
	 * Resolves MSFragger, FlashLFQ, optional Crystal-C / MSBooster / Sage /
	 * msconvert, and the marker FASTA. Mutates {@code params} with the path
	 * fields the headless runner reads.
	 */
	public static Tools prepare(DDAMagOpenSearchParameters params, MagDbItem mag, MetaResourceParameters mrp)
			throws IOException {
		return prepare(params, mag, mrp, true);
	}

	/** Evaluation runs stop at identification and do not require a FlashLFQ installation. */
	public static Tools prepareIdentification(DDAMagOpenSearchParameters params, MagDbItem mag,
			MetaResourceParameters mrp) throws IOException {
		return prepare(params, mag, mrp, false);
	}

	private static Tools prepare(DDAMagOpenSearchParameters params, MagDbItem mag, MetaResourceParameters mrp,
			boolean requireQuantification) throws IOException {
		if (params == null) {
			throw new IllegalArgumentException("params is null");
		}
		if (mag == null) {
			throw new IOException("MAG database item is required");
		}
		if (mrp == null) {
			throw new IOException("Resource parameters are required to locate FragPipe / FlashLFQ");
		}

		File msfraggerJar = mrp.resolveMsfraggerJar();
		if (msfraggerJar == null || !msfraggerJar.exists()) {
			throw new IOException("MSFragger.jar not found inside the configured FragPipe installation. "
					+ "Set FRAGPIPE to its bin/fragpipe launcher so the bundled jar can be detected.");
		}

		if (params.isStage1RescoringEnabled() && params.isStage1MsboosterEnabled()) {
			File msbooster = mrp.resolveMsboosterJar();
			File unimod = mrp.resolveUnimodObo();
			File diann = mrp.resolveMsboosterDiann();
			if (msbooster == null || diann == null || !diann.isFile()) {
				throw new IOException("Stage-1 MSBooster rescoring is enabled but "
						+ (msbooster == null
								? "MSBooster.jar was not found in the configured FragPipe installation's tools folder"
								: "no usable DIA-NN 1.8.x binary was found under tools/diann/. "
										+ "Expected diann, diann.exe, diann-linux*, or a versioned diann-* "
										+ "(e.g. diann-1.8.1.8) under tools/diann/<version>/ or "
										+ "tools/diann/<version>/<os>/. "
										+ diagnoseDiannLookup(mrp)
										+ "MSBooster drives DIA-NN through the 1.8.x prediction interface; "
										+ "a 2.x build silently predicts nothing")
						+ ". Set FRAGPIPE to the installation that bundles both, or switch stage-1 MSBooster off.");
			}
			params.setMsboosterJarPath(msbooster.getAbsolutePath());
			params.setDiannExecutablePath(diann.getAbsolutePath());
			if (unimod != null) {
				params.setUnimodOboPath(unimod.getAbsolutePath());
			}
		}

		if (params.isCrystalcEnabled()) {
			File crystalc = mrp.resolveCrystalcJar();
			if (crystalc == null) {
				throw new IOException("Crystal-C precursor refinement is enabled but no "
						+ "original-crystalc-*.jar was found in the configured FragPipe installation's "
						+ "tools folder. Set FRAGPIPE to an installation that bundles it, or switch the step off.");
			}
			params.setCrystalcJarPath(crystalc.getAbsolutePath());
		}

		String msconvertPath = mrp.getToolPath(ToolType.MSCONVERT);
		if (msconvertPath != null && !msconvertPath.isBlank()) {
			params.setMsconvertExecutable(msconvertPath);
		}

		File markerFasta = mag.getHapFastaFile(params.getGenomeSelectorMode());
		if (markerFasta == null || !markerFasta.exists()) {
			throw new IOException("MAG marker FASTA not found for genome selector mode "
					+ params.getGenomeSelectorMode());
		}

		String flashLfqPath = mrp.getToolPath(ToolType.FLASHLFQ);
		File flashLfqExe = flashLfqPath == null || flashLfqPath.isBlank() ? null : new File(flashLfqPath);
		if (requireQuantification && (flashLfqExe == null || !flashLfqExe.exists())) {
			throw new IOException("FlashLFQ not configured or missing. Set FLASHLFQ in the resources config "
					+ "to CMD.exe (Windows), CMD.dll (Linux/dotnet), or the FlashLFQ install directory.");
		}

		String funcDbPath = mrp.getToolPath(ToolType.FUNCDB);
		File funcDbFile = funcDbPath == null || funcDbPath.isBlank() ? null : new File(funcDbPath);
		if (funcDbFile != null && !funcDbFile.isFile()) {
			funcDbFile = null;
		}

		return new Tools(msfraggerJar, markerFasta, flashLfqExe, funcDbFile);
	}

	/**
	 * Short diagnostic for the MSBooster DIA-NN lookup failure: whether
	 * {@code tools/diann} exists and what version / OS folders were seen (so an
	 * empty {@code linux/} without a {@code diann} binary is obvious).
	 */
	static String diagnoseDiannLookup(MetaResourceParameters mrp) {
		if (mrp == null) {
			return "";
		}
		File tools = mrp.resolveFragpipeToolsDir();
		if (tools == null) {
			return "FragPipe tools/ was not resolved from FRAGPIPE. ";
		}
		File diannRoot = new File(tools, "diann");
		if (!diannRoot.isDirectory()) {
			return "Looked under " + tools.getAbsolutePath() + "/diann (missing). ";
		}
		File[] versions = diannRoot.listFiles(File::isDirectory);
		if (versions == null || versions.length == 0) {
			return "Found empty " + diannRoot.getAbsolutePath() + ". ";
		}
		StringBuilder sb = new StringBuilder("Looked under " + diannRoot.getAbsolutePath() + ": ");
		for (int i = 0; i < versions.length; i++) {
			if (i > 0) {
				sb.append("; ");
			}
			File v = versions[i];
			sb.append(v.getName());
			File[] kids = v.listFiles();
			if (kids == null || kids.length == 0) {
				sb.append(" (empty)");
				continue;
			}
			sb.append(" [");
			for (int j = 0; j < kids.length; j++) {
				if (j > 0) {
					sb.append(", ");
				}
				File k = kids[j];
				sb.append(k.getName());
				if (k.isDirectory()) {
					File[] nested = k.listFiles();
					int n = nested == null ? 0 : nested.length;
					sb.append('/').append(n).append(" entries");
				}
			}
			sb.append(']');
		}
		sb.append(". ");
		return sb.toString();
	}
}
