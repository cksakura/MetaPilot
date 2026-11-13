package uk.ac.quadram.metapilot.core.opensearch;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.function.Consumer;

import uk.ac.quadram.metapilot.core.variant.DiannFragmentPredictor;
import uk.ac.quadram.metapilot.core.variant.MzmlScanReader;

/**
 * Predicted-spectrum and retention-time agreement for the joint identification's candidates, from
 * DIA-NN 1.8.x predictions (the interface MSBooster and the variant form ranking already use) and
 * the run's mzML.
 * <p>
 * Per candidate:
 * <ul>
 * <li><b>ms2_entropy</b> - unweighted spectral-entropy similarity (Li et al. 2021; MSBooster's
 * {@code unweighted_spectral_entropy}) between DIA-NN's predicted b/y intensities and the most
 * intense observed peak within tolerance of each predicted fragment.</li>
 * <li><b>ms2_matched_fraction</b> - share of predicted intensity with an observed peak.</li>
 * <li><b>rt_error_norm</b> - |observed RT - calibrated predicted RT| in robust standard deviations,
 * calibrated per file.</li>
 * </ul>
 * The prediction is always of the unmodified backbone (C carbamidomethyl), so no modification needs
 * a DIA-NN name; masses are placed on the fragments afterwards. A closed candidate's variable
 * modifications sit where the closed search put them. An open candidate's delta sits on MSFragger's
 * best-localized residue, or, when it is not localized, on the residue that best explains the
 * spectrum. That choice is made the same way for targets and decoys, so it biases neither.
 * <p>
 * Candidates outside the closed search and beyond open rank {@code maxOpenRank} keep NaN: they are
 * rarely spectrum winners, and predicting them would multiply the cost.
 */
public final class PredictedSpectrumFeatures {

	static final double PROTON = 1.007276466812;
	static final double H2O = 18.010564684;
	private static final double LN4 = Math.log(4.0);
	/** Below this |delta| (Da) an open candidate carries no shift on its fragments. */
	static final double ZERO_DELTA_DA = 0.05;
	/** Most intense predicted fragments compared per candidate. */
	static final int TOP_FRAGMENTS = 24;
	/** Precursors per DIA-NN call; bounds the predictions held at once. */
	static final int CHUNK = 100_000;
	static final int MAX_RT_ANCHORS = 3000;

	public static final class Inputs {
		public File diann;
		public Map<String, File> mzmlBySample = new LinkedHashMap<>();
		public int threads = 1;
		public double fragmentTolerancePpm = 20.0;
		public int maxOpenRank = 2;
		public boolean carbamidomethyl = true;
		public File workDir;
	}

	private PredictedSpectrumFeatures() {
	}

	/**
	 * @return the per-candidate values; null when DIA-NN or the spectra are unavailable, with the
	 *         reason reported through {@code warn}
	 */
	public static JointFeatures.SpectrumScores compute(List<JointCandidate> candidates, Inputs in,
			Consumer<String> info, Consumer<String> warn) {
		if (in == null || in.diann == null || !in.diann.isFile()) {
			warn.accept("Joint identification: no DIA-NN 1.8.x executable configured; the predicted-spectrum and "
					+ "retention-time features are not used");
			return null;
		}
		Map<String, List<Integer>> byFile = new TreeMap<>();
		for (int i = 0; i < candidates.size(); i++) {
			JointCandidate c = candidates.get(i);
			if (eligible(c, in.maxOpenRank)) {
				byFile.computeIfAbsent(c.getFileName(), k -> new ArrayList<>()).add(i);
			}
		}
		Set<String> noSpectra = new HashSet<>(byFile.keySet());
		noSpectra.removeAll(in.mzmlBySample.keySet());
		if (noSpectra.size() == byFile.size()) {
			warn.accept("Joint identification: no mzML found for any run " + byFile.keySet()
					+ "; the predicted-spectrum features need <sample>.mzML beside the input spectra");
			return null;
		}
		if (!noSpectra.isEmpty()) {
			warn.accept("Joint identification: no mzML for " + noSpectra
					+ "; their candidates have no predicted-spectrum features");
		}
		JointFeatures.SpectrumScores out = new JointFeatures.SpectrumScores(candidates.size());
		int scored = 0;
		for (Map.Entry<String, List<Integer>> e : byFile.entrySet()) {
			File mzml = in.mzmlBySample.get(e.getKey());
			if (mzml == null) {
				continue;
			}
			try {
				scored += scoreFile(candidates, e.getValue(), mzml, in, out, info, warn);
			} catch (IOException | RuntimeException ex) {
				warn.accept("Joint identification: predicted-spectrum features failed for " + e.getKey()
						+ ", its candidates have none: " + ex.getMessage());
			}
		}
		int withRt = 0;
		for (float v : out.rtErrorNorm) {
			if (!Float.isNaN(v)) {
				withRt++;
			}
		}
		info.accept(String.format(Locale.ROOT,
				"Joint identification: predicted-spectrum features for %,d of %,d candidates, retention-time "
						+ "error for %,d",
				scored, candidates.size(), withRt));
		return scored == 0 && withRt == 0 ? null : out;
	}

	/** The cascade classifiers' names for the values {@link #attach} writes as external features. */
	public static final String ENTROPY = "ms2_entropy";
	public static final String MATCHED_FRACTION = "ms2_matched_fraction";
	public static final String ENTROPY_MARGIN = "ms2_entropy_margin";
	public static final String RT_ERROR = "rt_error_norm";
	public static final String[] PSM_FEATURES = { ENTROPY, MATCHED_FRACTION, ENTROPY_MARGIN, RT_ERROR };

	/**
	 * Predicted-spectrum features for one search's PSMs (the closed search's, or the open search's),
	 * written onto each PSM's external features under {@link #PSM_FEATURES}. A closed PSM is matched
	 * with its variable modifications, an open one with its delta placed as in {@link #score}.
	 * {@code ms2_entropy_margin} is the PSM's entropy minus the best other entropy of the same
	 * spectrum in this list (absent when the spectrum has no other scored candidate), so re-ranking
	 * among a spectrum's candidates can use it. Values that cannot be computed stay absent (missing).
	 *
	 * @return PSMs given an entropy; 0 when DIA-NN or the spectra are unavailable
	 */
	public static int attach(List<OpenSearchPsm> psms, Inputs in, String label, Consumer<String> info,
			Consumer<String> warn) {
		List<JointCandidate> candidates = JointCandidate.merge(psms, null, new int[2]);
		JointFeatures.SpectrumScores s = compute(candidates, in, m -> info.accept(label + ": " + m),
				m -> warn.accept(label + ": " + m));
		if (s == null) {
			return 0;
		}
		Map<String, double[]> bestTwo = new java.util.HashMap<>();
		for (int i = 0; i < candidates.size(); i++) {
			float e = s.entropy[i];
			if (Float.isNaN(e)) {
				continue;
			}
			double[] b = bestTwo.computeIfAbsent(candidates.get(i).spectrumKey(),
					k -> new double[] { Double.NEGATIVE_INFINITY, Double.NEGATIVE_INFINITY });
			if (e > b[0]) {
				b[1] = b[0];
				b[0] = e;
			} else if (e > b[1]) {
				b[1] = e;
			}
		}
		int withEntropy = 0;
		for (int i = 0; i < candidates.size(); i++) {
			JointCandidate c = candidates.get(i);
			Map<String, Double> values = new java.util.HashMap<>();
			float e = s.entropy[i];
			if (!Float.isNaN(e)) {
				values.put(ENTROPY, (double) e);
				values.put(MATCHED_FRACTION, (double) s.matchedFraction[i]);
				double[] b = bestTwo.get(c.spectrumKey());
				double other = e == b[0] ? b[1] : b[0];
				if (Double.isFinite(other)) {
					values.put(ENTROPY_MARGIN, e - other);
				}
				withEntropy++;
			}
			if (!Float.isNaN(s.rtErrorNorm[i])) {
				values.put(RT_ERROR, (double) s.rtErrorNorm[i]);
			}
			if (values.isEmpty()) {
				continue;
			}
			for (OpenSearchPsm psm : new OpenSearchPsm[] { c.closed, c.fragger, c.sage }) {
				if (psm != null) {
					Map<String, Double> external = new java.util.HashMap<>(psm.getExternalFeatures());
					external.putAll(values);
					psm.setExternalFeatures(external);
				}
			}
		}
		return withEntropy;
	}

	static boolean eligible(JointCandidate c, int maxOpenRank) {
		if (!standard(c.getPeptide()) || c.representative().getCharge() <= 0) {
			return false;
		}
		return c.closed != null || (c.fragger != null && c.fragger.getHitRank() <= maxOpenRank)
				|| (c.sage != null && c.sage.getHitRank() <= maxOpenRank);
	}

	private static int scoreFile(List<JointCandidate> candidates, List<Integer> rows, File mzml, Inputs in,
			JointFeatures.SpectrumScores out, Consumer<String> info, Consumer<String> warn) throws IOException {
		Set<Integer> scans = new HashSet<>();
		for (int i : rows) {
			scans.add(candidates.get(i).getScan());
		}
		info.accept(String.format(Locale.ROOT, "Joint identification: reading %,d MS2 scans from %s", scans.size(),
				mzml.getName()));
		MzmlScanReader.Result spectra = MzmlScanReader.read(mzml, scans);
		if (spectra.unsupportedEncoding > 0) {
			warn.accept("Joint identification: " + mzml.getName() + ": " + spectra.unsupportedEncoding
					+ " spectra use MS-Numpress and get no predicted-spectrum features");
		}
		Map<String, List<Integer>> byPrecursor = new LinkedHashMap<>();
		for (int i : rows) {
			JointCandidate c = candidates.get(i);
			byPrecursor.computeIfAbsent(DiannFragmentPredictor.key(c.getPeptide(), c.representative().getCharge()),
					k -> new ArrayList<>()).add(i);
		}
		List<String> keys = new ArrayList<>(byPrecursor.keySet());
		double[] iRT = new double[rows.size()];
		Arrays.fill(iRT, Double.NaN);
		Map<Integer, Integer> slot = new java.util.HashMap<>(rows.size() * 2);
		for (int k = 0; k < rows.size(); k++) {
			slot.put(rows.get(k), k);
		}
		File work = new File(in.workDir, safe(mzml.getName()));
		int scored = 0;
		for (int from = 0; from < keys.size(); from += CHUNK) {
			List<String> chunk = keys.subList(from, Math.min(keys.size(), from + CHUNK));
			List<String> peptides = new ArrayList<>(chunk.size());
			List<Integer> charges = new ArrayList<>(chunk.size());
			for (String key : chunk) {
				int cut = key.lastIndexOf('/');
				peptides.add(key.substring(0, cut));
				charges.add(Integer.parseInt(key.substring(cut + 1)));
			}
			Map<String, DiannFragmentPredictor.Prediction> predicted;
			try {
				predicted = DiannFragmentPredictor.predict(in.diann, peptides, charges, in.carbamidomethyl, work,
						in.threads, info);
			} catch (IOException ex) {
				// One chunk DIA-NN cannot predict (it returns fewer records than it was given, say) costs
				// that chunk only; the file's other chunks and its retention-time calibration go on.
				warn.accept(String.format(Locale.ROOT,
						"Joint identification: DIA-NN failed for precursors %,d-%,d of %s; they have no "
								+ "predicted-spectrum features: %s",
						from + 1, from + chunk.size(), mzml.getName(), ex.getMessage()));
				continue;
			}
			for (String key : chunk) {
				DiannFragmentPredictor.Prediction p = predicted.get(key);
				if (p == null) {
					continue;
				}
				for (int i : byPrecursor.get(key)) {
					JointCandidate c = candidates.get(i);
					iRT[slot.get(i)] = p.iRT;
					MzmlScanReader.Peaks peaks = spectra.scans.get(c.getScan());
					if (peaks == null) {
						continue;
					}
					double[] s = score(c, p, peaks, in.fragmentTolerancePpm);
					if (!Double.isNaN(s[0])) {
						out.entropy[i] = (float) s[0];
						out.matchedFraction[i] = (float) s[1];
						scored++;
					}
				}
			}
		}
		calibrateRetentionTime(candidates, rows, iRT, out);
		return scored;
	}

	private static String safe(String name) {
		return name.replaceAll("[^A-Za-z0-9._-]", "_");
	}

	// ------------------------------------------------------------ MS2 similarity

	/** {entropy similarity, matched predicted-intensity fraction}; NaN when nothing could be compared. */
	static double[] score(JointCandidate c, DiannFragmentPredictor.Prediction p, MzmlScanReader.Peaks peaks,
			double tolerancePpm) {
		String peptide = c.getPeptide();
		int n = peptide.length();
		int charge = c.representative().getCharge();
		double[] mods = new double[n + 1];
		int[] candidateSites;
		double delta = 0.0;
		if (c.closed != null) {
			double[] closedMods = c.closed.getClosedModifications();
			if (closedMods != null) {
				for (int k = 0; k < Math.min(closedMods.length, mods.length); k++) {
					mods[k] = closedMods[k];
				}
			}
			candidateSites = new int[] { -1 };
		} else {
			OpenSearchPsm open = c.openPrimary();
			delta = open.getNormalizedMassDiff();
			if (!Double.isFinite(delta) || Math.abs(delta) < ZERO_DELTA_DA) {
				delta = 0.0;
				candidateSites = new int[] { -1 };
			} else {
				int site = open.getFirstLocalizedSiteIndex();
				if (site >= 0 && site < n) {
					candidateSites = new int[] { site };
				} else {
					candidateSites = new int[n];
					for (int k = 0; k < n; k++) {
						candidateSites[k] = k;
					}
				}
			}
		}
		int[] top = topFragments(p, n, charge);
		if (top.length == 0) {
			return new double[] { Double.NaN, Double.NaN };
		}
		double bestEntropy = Double.NaN;
		double bestMatched = Double.NaN;
		double[] predicted = new double[top.length];
		double[] observed = new double[top.length];
		for (int site : candidateSites) {
			if (site >= 0) {
				Arrays.fill(mods, 0.0);
				mods[site + 1] = delta;
			}
			double matched = 0.0;
			double total = 0.0;
			for (int k = 0; k < top.length; k++) {
				int f = top[k];
				double mz = fragmentMz(peptide, mods, p.ionNumber[f], p.yIon[f], p.charge[f]);
				double tol = mz * tolerancePpm * 1e-6;
				predicted[k] = p.intensity[f];
				observed[k] = peaks.maxIntensity(mz - tol, mz + tol);
				total += predicted[k];
				if (observed[k] > 0) {
					matched += predicted[k];
				}
			}
			double e = entropySimilarity(predicted, observed, top.length);
			if (Double.isNaN(bestEntropy) || e > bestEntropy) {
				bestEntropy = e;
				bestMatched = total > 0 ? matched / total : 0.0;
			}
		}
		return new double[] { bestEntropy, bestMatched };
	}

	/** Indices of the most intense usable predicted fragments, most intense first. */
	static int[] topFragments(DiannFragmentPredictor.Prediction p, int length, int precursorCharge) {
		List<Integer> usable = new ArrayList<>();
		for (int f = 0; f < p.size(); f++) {
			if (p.intensity[f] > 0 && p.ionNumber[f] >= 1 && p.ionNumber[f] < length
					&& p.charge[f] <= Math.max(1, precursorCharge)) {
				usable.add(f);
			}
		}
		usable.sort((a, b) -> Float.compare(p.intensity[b], p.intensity[a]));
		int m = Math.min(TOP_FRAGMENTS, usable.size());
		int[] out = new int[m];
		for (int k = 0; k < m; k++) {
			out[k] = usable.get(k);
		}
		return out;
	}

	/**
	 * Fragment m/z with per-position modification masses ({@code mods[0]} the N-terminus,
	 * {@code mods[i]} residue {@code i}, 1-based). b ions hold residues 1..ion and the N-terminus, y
	 * ions the last {@code ion} residues.
	 */
	static double fragmentMz(String peptide, double[] mods, int ion, boolean yIon, int charge) {
		int n = peptide.length();
		int from = yIon ? n - ion : 0;
		int to = yIon ? n : ion;
		double mass = yIon ? H2O : mods[0];
		for (int i = from; i < to; i++) {
			mass += SubstitutionMassTable.residueMass(peptide.charAt(i)) + mods[i + 1];
		}
		return (mass + charge * PROTON) / charge;
	}

	/** {@code 1 - (2 H(P+Q) - H(P) - H(Q)) / ln 4}, each vector normalised to 1; 0 when either is empty. */
	public static double entropySimilarity(double[] p, double[] q, int n) {
		double sp = 0;
		double sq = 0;
		for (int i = 0; i < n; i++) {
			sp += p[i];
			sq += q[i];
		}
		if (sp <= 0 || sq <= 0) {
			return 0.0;
		}
		double hp = 0;
		double hq = 0;
		double hm = 0;
		for (int i = 0; i < n; i++) {
			double a = p[i] / sp;
			double b = q[i] / sq;
			double m = (a + b) / 2;
			hp -= xlogx(a);
			hq -= xlogx(b);
			hm -= xlogx(m);
		}
		return Math.max(0.0, Math.min(1.0, 1.0 - (2 * hm - hp - hq) / LN4));
	}

	private static double xlogx(double x) {
		return x > 0 ? x * Math.log(x) : 0.0;
	}

	static boolean standard(String peptide) {
		if (peptide == null || peptide.isEmpty()) {
			return false;
		}
		for (int i = 0; i < peptide.length(); i++) {
			if ("ACDEFGHIKLMNPQRSTVWY".indexOf(peptide.charAt(i)) < 0) {
				return false;
			}
		}
		return true;
	}

	// ------------------------------------------------------------ retention time

	/**
	 * Per file, maps DIA-NN's iRT onto the observed retention time with a monotone piecewise-linear
	 * curve through binned medians of anchor candidates, and writes each candidate's absolute error
	 * in robust standard deviations of the anchors' residuals.
	 * <p>
	 * Anchors are the file's candidates with the strongest search evidence (closed or open rank 1,
	 * best {@code -log10(expect)}), chosen by score alone - never by target/decoy label - so the
	 * calibration treats a decoy exactly as it treats a target. A wrong anchor moves a bin median
	 * only if it is a large share of its bin.
	 */
	static void calibrateRetentionTime(List<JointCandidate> candidates, List<Integer> rows, double[] iRT,
			JointFeatures.SpectrumScores out) {
		List<Integer> anchors = new ArrayList<>();
		for (int k = 0; k < rows.size(); k++) {
			JointCandidate c = candidates.get(rows.get(k));
			if (Double.isFinite(iRT[k]) && rankOne(c) && Double.isFinite(c.representative().getRetentionTime())) {
				anchors.add(k);
			}
		}
		if (anchors.size() < 50) {
			return;
		}
		anchors.sort((a, b) -> Double.compare(evidence(candidates.get(rows.get(b))),
				evidence(candidates.get(rows.get(a)))));
		List<Integer> top = anchors.subList(0, Math.min(MAX_RT_ANCHORS, anchors.size()));
		double[][] points = new double[top.size()][2];
		for (int k = 0; k < top.size(); k++) {
			points[k][0] = iRT[top.get(k)];
			points[k][1] = candidates.get(rows.get(top.get(k))).representative().getRetentionTime();
		}
		Arrays.sort(points, (a, b) -> Double.compare(a[0], b[0]));
		int bins = Math.max(4, Math.min(30, points.length / 50));
		double[] bx = new double[bins];
		double[] by = new double[bins];
		for (int b = 0; b < bins; b++) {
			int from = b * points.length / bins;
			int to = (b + 1) * points.length / bins;
			double[] xs = new double[to - from];
			double[] ys = new double[to - from];
			for (int k = from; k < to; k++) {
				xs[k - from] = points[k][0];
				ys[k - from] = points[k][1];
			}
			bx[b] = median(xs);
			by[b] = median(ys);
		}
		isotonic(by);
		double[] residuals = new double[points.length];
		for (int k = 0; k < points.length; k++) {
			residuals[k] = Math.abs(points[k][1] - interpolate(bx, by, points[k][0]));
		}
		double scale = 1.4826 * median(residuals);
		if (!(scale > 0)) {
			return;
		}
		for (int k = 0; k < rows.size(); k++) {
			if (!Double.isFinite(iRT[k])) {
				continue;
			}
			double rt = candidates.get(rows.get(k)).representative().getRetentionTime();
			if (Double.isFinite(rt)) {
				out.rtErrorNorm[rows.get(k)] = (float) (Math.abs(rt - interpolate(bx, by, iRT[k])) / scale);
			}
		}
	}

	private static boolean rankOne(JointCandidate c) {
		return (c.closed != null && c.closed.getHitRank() == 1) || (c.fragger != null && c.fragger.getHitRank() == 1);
	}

	private static double evidence(JointCandidate c) {
		double best = Double.NEGATIVE_INFINITY;
		for (OpenSearchPsm psm : new OpenSearchPsm[] { c.closed, c.fragger }) {
			if (psm != null) {
				double v = JointFeatures.negLog10(psm.getExpect());
				if (!Double.isNaN(v)) {
					best = Math.max(best, v);
				}
			}
		}
		return best;
	}

	/** Pool-adjacent-violators, in place: the fitted curve may not run backwards in time. */
	public static void isotonic(double[] y) {
		int n = y.length;
		double[] value = new double[n];
		int[] weight = new int[n];
		int blocks = 0;
		for (int i = 0; i < n; i++) {
			value[blocks] = y[i];
			weight[blocks] = 1;
			blocks++;
			while (blocks > 1 && value[blocks - 2] > value[blocks - 1]) {
				double merged = (value[blocks - 2] * weight[blocks - 2] + value[blocks - 1] * weight[blocks - 1])
						/ (weight[blocks - 2] + weight[blocks - 1]);
				weight[blocks - 2] += weight[blocks - 1];
				value[blocks - 2] = merged;
				blocks--;
			}
		}
		int i = 0;
		for (int b = 0; b < blocks; b++) {
			for (int k = 0; k < weight[b]; k++) {
				y[i++] = value[b];
			}
		}
	}

	/** Linear interpolation through (x, y), extended linearly beyond the ends. */
	public static double interpolate(double[] x, double[] y, double v) {
		int n = x.length;
		if (n == 1) {
			return y[0];
		}
		int hi = 1;
		while (hi < n - 1 && x[hi] < v) {
			hi++;
		}
		int lo = hi - 1;
		double dx = x[hi] - x[lo];
		if (dx <= 0) {
			return (y[lo] + y[hi]) / 2;
		}
		return y[lo] + (v - x[lo]) * (y[hi] - y[lo]) / dx;
	}

	public static double median(double[] values) {
		double[] v = values.clone();
		Arrays.sort(v);
		int n = v.length;
		if (n == 0) {
			return Double.NaN;
		}
		return n % 2 == 1 ? v[n / 2] : (v[n / 2 - 1] + v[n / 2]) / 2;
	}
}
