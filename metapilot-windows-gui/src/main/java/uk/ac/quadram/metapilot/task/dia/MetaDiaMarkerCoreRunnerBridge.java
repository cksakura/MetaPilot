package uk.ac.quadram.metapilot.task.dia;

import java.io.File;
import java.util.Arrays;

import org.apache.logging.log4j.Logger;

import uk.ac.quadram.metapilot.core.report.PythonEnvironmentResolver;
import uk.ac.quadram.metapilot.core.workflow.ProgressCallback;
import uk.ac.quadram.metapilot.resources.MetaResourceParameters;
import uk.ac.quadram.metapilot.resources.MetaResourceParameters.ToolType;
import uk.ac.quadram.metapilot.task.dia.par.MetaParameterDiaMarker;
import uk.ac.quadram.metapilot.task.par.MetaData;

final class MetaDiaMarkerCoreRunnerBridge {

	@FunctionalInterface
	interface ProgressAction {
		void update(double workDone, double max);
	}

	@FunctionalInterface
	interface MessageAction {
		void update(String message);
	}

	private MetaDiaMarkerCoreRunnerBridge() {
	}

	static ProgressCallback createGuiProgressCallback(ProgressAction progressAction,
			MessageAction messageAction, Logger logger, String taskName) {
		return new ProgressCallback() {
			@Override
			public void onProgress(double workDone, double max) {
				if (max > 0) {
					progressAction.update(workDone, max);
				}
			}

			@Override
			public void onMessage(String message) {
				if (message != null && !message.isBlank()) {
					messageAction.update(message);
				}
			}

			@Override
			public void onComplete(Object result) {
			}

			@Override
			public void onError(Exception error) {
				if (error != null && logger != null) {
					logger.error(taskName + ": core runner error", error);
				}
			}
		};
	}

	static MetaParameterDiaMarker buildCoreRunnerParams(MetaParameterDiaMarker original,
			MetaData metadata, String[] diaFiles) {
		MetaParameterDiaMarker copy = new MetaParameterDiaMarker(
				original,
				original.getMagDbItem(),
				original.isLibrarySearch(),
				original.isSelfModelSearch(),
				original.getLibrary(),
				original.isUseModel(),
				original.getModels());
		copy.setMarkerMode(original.getMarkerMode());
		copy.setGenomesPerRound(original.getGenomesPerRound());
		copy.setMExcision(original.isMExcision());
		copy.setCarbamiC(original.isCarbamiC());
		copy.setOxidationM(original.isOxidationM());
		copy.setAcetyProN(original.isAcetyProN());
		copy.setPhosphoSTY(original.isPhosphoSTY());
		copy.setKGG(original.isKGG());
		copy.setUsePeptideBasedRound(original.isUsePeptideBasedRound());
		copy.setPeptideDetectabilityThreshold(original.getPeptideDetectabilityThreshold());
		copy.setMetadata(cloneMetadataWithRawFiles(metadata, diaFiles));
		copy.setDbSearchResultFile(original.getDbSearchResultFile());
		copy.setReportDir(original.getReportDir());
		return copy;
	}

	static File resolveGuiToolFile(MetaResourceParameters mrp, ToolType toolType) {
		if (mrp == null || toolType == null) {
			return null;
		}
		String path = mrp.getToolPath(toolType);
		if (path == null || path.isBlank()) {
			return null;
		}
		return new File(path);
	}

	static File inferCondaPrefixFromPython(File pythonExe) {
		return PythonEnvironmentResolver.inferCondaPrefixFromPython(pythonExe);
	}

	private static MetaData cloneMetadataWithRawFiles(MetaData source, String[] rawFilesForRunner) {
		MetaData md = new MetaData();
		if (source == null) {
			md.setRawFiles(rawFilesForRunner == null ? new String[0] : rawFilesForRunner.clone());
			return md;
		}
		md.setRawFiles(rawFilesForRunner == null ? new String[0] : rawFilesForRunner.clone());
		md.setExpNames(clone(source.getExpNames()));
		md.setFractions(clone(source.getFractions()));
		md.setReplicates(clone(source.getReplicates()));
		md.setMetaTypeCount(source.getMetaTypeCount());
		md.setMetaInfo(clone2d(source.getMetaInfo()));
		md.setMetaTypeTitles(clone(source.getMetaTypeTitles()));
		md.setLabelTitle(clone(source.getLabelTitle()));
		md.setLabelExpNames(clone(source.getLabelExpNames()));
		md.setIsobaricReference(clone(source.getIsobaricReference()));
		md.setSelectRef(clone(source.getSelectRef()));
		md.setRefChannelId(source.getRefChannelId());
		return md;
	}

	private static String[] clone(String[] values) {
		return values == null ? null : Arrays.copyOf(values, values.length);
	}

	private static int[] clone(int[] values) {
		return values == null ? null : Arrays.copyOf(values, values.length);
	}

	private static boolean[] clone(boolean[] values) {
		return values == null ? null : Arrays.copyOf(values, values.length);
	}

	private static String[][] clone2d(String[][] values) {
		if (values == null) {
			return null;
		}
		String[][] copy = new String[values.length][];
		for (int i = 0; i < values.length; i++) {
			copy[i] = clone(values[i]);
		}
		return copy;
	}
}
