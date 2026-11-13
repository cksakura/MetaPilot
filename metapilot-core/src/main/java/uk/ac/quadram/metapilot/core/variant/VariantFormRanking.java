package uk.ac.quadram.metapilot.core.variant;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.EnumMap;
import java.util.HashSet;
import java.util.IdentityHashMap;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.function.Consumer;

/**
 * MS2 form ranking for open-search candidates: does the spectrum prefer the
 * candidate's substitution over the same delta mass as a modification or as a
 * different substitution? Sets each open-search candidate's
 * {@link VariantCandidate#ms2Verdict} and scores; pFind candidates are left at
 * NOT_APPLICABLE, their calls trusted as-is.
 * <p>
 * Report-only and secondary: the verdict and scores are columns of
 * {@code variant_candidates.tsv}; they do not set the candidate's class or the
 * summary. Never throws; whatever cannot
 * be ranked (ranking off, no DIA-NN, no mzML, scan missing, prediction failed)
 * stays NOT_RANKED with the reason logged.
 */
public final class VariantFormRanking {

	private static final String LOG_PREFIX = "Pangenome variant check, MS2 form ranking: ";
	static final String WORK_SUBDIR = "variant_ranking";

	/** What the ranking needs from the workflow. */
	public static final class Inputs {
		public boolean enabled = true;
		/** DIA-NN 1.8.x executable; ranking is skipped without it. */
		public File diann;
		public int threads = 1;
		/** The run's spectrum files; a non-mzML entry is looked for as {@code <base>.mzML} beside it. */
		public List<File> spectraFiles = new ArrayList<>();
		public double margin = VariantFormRankingSettings.DEFAULT_MARGIN;
		public double fragmentTolerancePpm = VariantFormRankingSettings.DEFAULT_FRAGMENT_TOLERANCE_PPM;
		/** Tolerance for the same-mass substitutions competing with the candidate's own. */
		public double substitutionToleranceDa = 0.005;
		/** Whether C carries the fixed carbamidomethyl in the search. */
		public boolean carbamidomethyl = true;

		public void apply(VariantFormRankingSettings settings) {
			enabled = settings.isEnabled();
			margin = settings.getMargin();
			fragmentTolerancePpm = settings.getFragmentTolerancePpm();
		}
	}

	private VariantFormRanking() {
	}

	static void run(List<VariantCandidate> candidates, Inputs in, File outputDir, Consumer<String> info,
			Consumer<String> warn) {
		List<VariantCandidate> open = new ArrayList<>();
		for (VariantCandidate c : candidates) {
			if (c.engine == VariantCandidate.Engine.OPEN_SEARCH) {
				c.ms2Verdict = VariantCandidate.Ms2Verdict.NOT_RANKED;
				open.add(c);
			}
		}
		if (open.isEmpty()) {
			return;
		}
		if (in == null || !in.enabled) {
			info.accept(LOG_PREFIX + (in == null ? "not configured" : "off (" + VariantFormRankingSettings.KEY_ENABLED
					+ " = false)") + "; " + open.size() + " candidates NOT_RANKED");
			return;
		}
		if (in.diann == null || !in.diann.isFile()) {
			warn.accept(LOG_PREFIX + "no DIA-NN executable configured" + (in.diann == null ? "" : " (" + in.diann
					+ " not found)") + "; it predicts the fragment intensities that separate a substitution from "
					+ "a same-site modification. " + open.size() + " candidates NOT_RANKED");
			return;
		}
		try {
			rank(open, in, new File(outputDir, WORK_SUBDIR), info, warn);
		} catch (IOException | RuntimeException e) {
			for (VariantCandidate c : open) {
				clear(c);
			}
			warn.accept(LOG_PREFIX + "failed, candidates NOT_RANKED: " + e.getMessage());
		}
		Map<VariantCandidate.Ms2Verdict, Integer> counts = new EnumMap<>(VariantCandidate.Ms2Verdict.class);
		int elsewhere = 0;
		for (VariantCandidate c : open) {
			counts.merge(c.ms2Verdict, 1, Integer::sum);
			if (Boolean.TRUE.equals(c.ms2LocalisesElsewhere)) {
				elsewhere++;
			}
		}
		info.accept(LOG_PREFIX + "same-site verdicts " + counts + ", " + elsewhere
				+ " localise elsewhere (margin " + in.margin + ", provisional)");
	}

	private static void rank(List<VariantCandidate> open, Inputs in, File workDir, Consumer<String> info,
			Consumer<String> warn) throws IOException {
		Map<String, File> mzml = mzmlBySample(in.spectraFiles);
		Map<String, List<VariantCandidate>> byFile = new TreeMap<>();
		Map<VariantCandidate, List<VariantFormRanker.Form>> forms = new IdentityHashMap<>();
		Set<String> keys = new LinkedHashSet<>();
		List<String> peptides = new ArrayList<>();
		List<Integer> charges = new ArrayList<>();
		int noCharge = 0;
		int noForms = 0;
		Set<String> noSpectra = new HashSet<>();
		for (VariantCandidate c : open) {
			if (c.charge <= 0) {
				noCharge++;
				continue;
			}
			String sample = c.file == null ? "" : c.file;
			if (!mzml.containsKey(sample)) {
				noSpectra.add(sample);
				continue;
			}
			List<VariantFormRanker.Form> f = VariantFormRanker.forms(c, in.substitutionToleranceDa);
			if (f.isEmpty()) {
				noForms++;
				continue;
			}
			forms.put(c, f);
			byFile.computeIfAbsent(sample, k -> new ArrayList<>()).add(c);
			for (VariantFormRanker.Form form : f) {
				if (keys.add(DiannFragmentPredictor.key(form.backbone, c.charge))) {
					peptides.add(form.backbone);
					charges.add(c.charge);
				}
			}
		}
		if (noCharge > 0) {
			warn.accept(LOG_PREFIX + noCharge + " candidates without a precursor charge, NOT_RANKED");
		}
		if (noForms > 0) {
			warn.accept(LOG_PREFIX + noForms + " candidates with a non-standard residue, NOT_RANKED");
		}
		if (!noSpectra.isEmpty()) {
			warn.accept(LOG_PREFIX + "no mzML found for " + noSpectra + "; their candidates are NOT_RANKED. "
					+ "Ranking reads the MS2 peaks from <sample>.mzML beside the input spectra");
		}
		if (peptides.isEmpty()) {
			return;
		}
		info.accept(String.format(Locale.ROOT, LOG_PREFIX + "%,d candidates in %d files, %,d forms to predict",
				forms.size(), byFile.size(), peptides.size()));

		Map<String, DiannFragmentPredictor.Prediction> predicted = DiannFragmentPredictor.predict(in.diann, peptides,
				charges, in.carbamidomethyl, workDir, in.threads, info);
		VariantFormRanker.Predictions lookup = (backbone, charge) -> predicted
				.get(DiannFragmentPredictor.key(backbone, charge));

		int missingScans = 0;
		for (Map.Entry<String, List<VariantCandidate>> e : byFile.entrySet()) {
			if (Thread.currentThread().isInterrupted()) {
				throw new IOException("interrupted");
			}
			Set<Integer> wanted = new HashSet<>();
			for (VariantCandidate c : e.getValue()) {
				wanted.add(c.scan);
			}
			File file = mzml.get(e.getKey());
			MzmlScanReader.Result spectra;
			try {
				spectra = MzmlScanReader.read(file, wanted);
			} catch (IOException ex) {
				warn.accept(LOG_PREFIX + "cannot read " + file + ", its candidates are NOT_RANKED: " + ex.getMessage());
				continue;
			}
			if (spectra.unsupportedEncoding > 0) {
				warn.accept(LOG_PREFIX + file.getName() + ": " + spectra.unsupportedEncoding
						+ " spectra use MS-Numpress, which is not decoded; their candidates are NOT_RANKED");
			}
			for (VariantCandidate c : e.getValue()) {
				MzmlScanReader.Peaks peaks = spectra.scans.get(c.scan);
				if (peaks == null) {
					missingScans++;
					continue;
				}
				VariantFormRanker.rank(c, forms.get(c), lookup, peaks, in.fragmentTolerancePpm, in.margin);
			}
		}
		if (missingScans > 0) {
			warn.accept(LOG_PREFIX + missingScans + " candidates whose scan was not found as an MS2 spectrum in "
					+ "its mzML, NOT_RANKED");
		}
	}

	private static void clear(VariantCandidate c) {
		c.ms2Verdict = VariantCandidate.Ms2Verdict.NOT_RANKED;
		c.ms2ScoreSubstitution = Double.NaN;
		c.ms2ScoreOtherSubstitution = Double.NaN;
		c.ms2ScoreModification = Double.NaN;
		c.ms2ScoreElsewhere = Double.NaN;
		c.ms2LocalisesElsewhere = null;
		c.ms2BestForm = "";
	}

	/** Sample base name to its mzML; open-search PSMs carry the base name of the spectrum file. */
	public static Map<String, File> mzmlBySample(List<File> spectraFiles) {
		Map<String, File> out = new LinkedHashMap<>();
		if (spectraFiles == null) {
			return out;
		}
		for (File f : spectraFiles) {
			if (f == null) {
				continue;
			}
			String name = f.getName();
			int dot = name.lastIndexOf('.');
			String base = dot > 0 ? name.substring(0, dot) : name;
			File candidate = name.toLowerCase(Locale.ROOT).endsWith(".mzml") ? f : findMzml(f.getParentFile(), base);
			if (candidate != null && candidate.isFile()) {
				out.putIfAbsent(base, candidate);
			}
		}
		return out;
	}

	private static File findMzml(File dir, String base) {
		if (dir == null) {
			return null;
		}
		for (String ext : new String[] { ".mzML", ".mzml", ".MZML" }) {
			File f = new File(dir, base + ext);
			if (f.isFile()) {
				return f;
			}
		}
		return null;
	}
}
