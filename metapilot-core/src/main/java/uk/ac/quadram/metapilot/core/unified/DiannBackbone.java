package uk.ac.quadram.metapilot.core.unified;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.function.Consumer;

import uk.ac.quadram.metapilot.core.opensearch.PredictedSpectrumFeatures;
import uk.ac.quadram.metapilot.core.variant.DiannFragmentPredictor;

/**
 * DIA-NN 1.8.x prediction of the unmodified backbone (C carbamidomethyl only), matched to the
 * observed spectrum <b>without</b> the hypothesis's mass difference or variable modifications: a
 * fragment carrying the shifted residue simply does not match. That costs every shifted hypothesis
 * alike, and FDR is decided within a batch, where every hypothesis carries the same kind of shift.
 * <p>
 * Features: spectral-entropy similarity and Pearson correlation over the most intense predicted
 * fragments, the share of predicted intensity observed, and the retention-time error after a
 * per-run isotonic calibration on label-blind anchors (the strongest candidates by annotation
 * score).
 */
final class DiannBackbone {

	static final int TOP_FRAGMENTS = 24;
	static final int CHUNK = 100_000;
	static final int MAX_RT_ANCHORS = 3000;

	private DiannBackbone() {
	}

	/** Scores the hypotheses of one run whose spectra are in {@code spectra} (by spectrum index). */
	static void scoreRun(List<Hypothesis> hyps, Map<Integer, SpectrumAnnotator.Spectrum> spectra, File diann,
			File workDir, int threads, double tolerancePpm, Consumer<String> log, Consumer<String> warn) {
		Map<String, List<Hypothesis>> byPrecursor = new LinkedHashMap<>();
		for (Hypothesis h : hyps) {
			if (h.neighbour || !standard(h.peptide) || h.charge <= 0) {
				continue;
			}
			byPrecursor.computeIfAbsent(DiannFragmentPredictor.key(h.peptide, h.charge), k -> new ArrayList<>()).add(h);
		}
		List<String> keys = new ArrayList<>(byPrecursor.keySet());
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
				predicted = DiannFragmentPredictor.predict(diann, peptides, charges, true, workDir, threads, log);
			} catch (IOException e) {
				warn.accept(String.format(Locale.ROOT, "DIA-NN failed for precursors %,d-%,d: %s", from + 1,
						from + chunk.size(), e.getMessage()));
				continue;
			}
			for (String key : chunk) {
				DiannFragmentPredictor.Prediction p = predicted.get(key);
				if (p == null) {
					continue;
				}
				for (Hypothesis h : byPrecursor.get(key)) {
					SpectrumAnnotator.Spectrum s = spectra.get(h.spectrum);
					h.predIrt = p.iRT;
					if (s != null) {
						match(h, p, s, tolerancePpm);
					}
				}
			}
		}
		calibrateRt(hyps);
	}

	static void match(Hypothesis h, DiannFragmentPredictor.Prediction p, SpectrumAnnotator.Spectrum s,
			double tolerancePpm) {
		int n = h.peptide.length();
		SpectrumAnnotator.Ladder ladder = new SpectrumAnnotator.Ladder(h.peptide, Math.max(2, h.charge + 1));
		List<Integer> usable = new ArrayList<>();
		for (int f = 0; f < p.size(); f++) {
			if (p.intensity[f] > 0 && p.ionNumber[f] >= 1 && p.ionNumber[f] < n && p.charge[f] <= Math.max(1, h.charge)) {
				usable.add(f);
			}
		}
		usable.sort((a, b) -> Float.compare(p.intensity[b], p.intensity[a]));
		int m = Math.min(TOP_FRAGMENTS, usable.size());
		if (m == 0) {
			return;
		}
		double[] pred = new double[m];
		double[] obs = new double[m];
		double total = 0;
		double seen = 0;
		for (int k = 0; k < m; k++) {
			int f = usable.get(k);
			double neutral = p.yIon[f] ? ladder.prefix[n] - ladder.prefix[n - p.ionNumber[f]] + SpectrumAnnotator.H2O
					: ladder.prefix[p.ionNumber[f]];
			double mz = (neutral + p.charge[f] * SpectrumAnnotator.PROTON) / p.charge[f];
			double tol = mz * tolerancePpm * 1e-6;
			pred[k] = p.intensity[f];
			obs[k] = maxIntensity(s, mz - tol, mz + tol);
			total += pred[k];
			if (obs[k] > 0) {
				seen += pred[k];
			}
		}
		h.predEntropy = (float) PredictedSpectrumFeatures.entropySimilarity(pred, obs, m);
		h.predPearson = (float) pearson(pred, obs);
		h.predObserved = (float) (total > 0 ? seen / total : 0);
	}

	private static double maxIntensity(SpectrumAnnotator.Spectrum s, double lo, double hi) {
		int i = Arrays.binarySearch(s.mz, lo);
		if (i < 0) {
			i = -i - 1;
		}
		double best = 0;
		for (; i < s.mz.length && s.mz[i] <= hi; i++) {
			best = Math.max(best, s.intensity[i]);
		}
		return best;
	}

	static double pearson(double[] a, double[] b) {
		int n = a.length;
		double ma = 0;
		double mb = 0;
		for (int i = 0; i < n; i++) {
			ma += a[i];
			mb += b[i];
		}
		ma /= n;
		mb /= n;
		double sab = 0;
		double saa = 0;
		double sbb = 0;
		for (int i = 0; i < n; i++) {
			sab += (a[i] - ma) * (b[i] - mb);
			saa += (a[i] - ma) * (a[i] - ma);
			sbb += (b[i] - mb) * (b[i] - mb);
		}
		return saa <= 0 || sbb <= 0 ? 0 : sab / Math.sqrt(saa * sbb);
	}

	/**
	 * Per run: observed RT against predicted iRT through binned medians of the anchors (closed or
	 * open rank-1 candidates with the highest annotation score, chosen by score alone), made
	 * monotone; the error is in robust SDs of the anchors' residuals.
	 */
	static void calibrateRt(List<Hypothesis> hyps) {
		List<Hypothesis> anchors = new ArrayList<>();
		for (Hypothesis h : hyps) {
			if (!h.neighbour && Double.isFinite(h.predIrt) && Double.isFinite(h.retentionTime) && h.annotation != null
					&& (h.foundByClosed || h.openRank == 1)) {
				anchors.add(h);
			}
		}
		if (anchors.size() < 50) {
			return;
		}
		anchors.sort((a, b) -> Float.compare(b.annotation.eHyperscore, a.annotation.eHyperscore));
		List<Hypothesis> top = anchors.subList(0, Math.min(MAX_RT_ANCHORS, anchors.size()));
		double[][] pts = new double[top.size()][2];
		for (int k = 0; k < top.size(); k++) {
			pts[k][0] = top.get(k).predIrt;
			pts[k][1] = top.get(k).retentionTime;
		}
		Arrays.sort(pts, (a, b) -> Double.compare(a[0], b[0]));
		int bins = Math.max(4, Math.min(30, pts.length / 50));
		double[] bx = new double[bins];
		double[] by = new double[bins];
		for (int b = 0; b < bins; b++) {
			int from = b * pts.length / bins;
			int to = (b + 1) * pts.length / bins;
			double[] xs = new double[to - from];
			double[] ys = new double[to - from];
			for (int k = from; k < to; k++) {
				xs[k - from] = pts[k][0];
				ys[k - from] = pts[k][1];
			}
			bx[b] = PredictedSpectrumFeatures.median(xs);
			by[b] = PredictedSpectrumFeatures.median(ys);
		}
		PredictedSpectrumFeatures.isotonic(by);
		double[] res = new double[pts.length];
		for (int k = 0; k < pts.length; k++) {
			res[k] = Math.abs(pts[k][1] - PredictedSpectrumFeatures.interpolate(bx, by, pts[k][0]));
		}
		double scale = 1.4826 * PredictedSpectrumFeatures.median(res);
		if (!(scale > 0)) {
			return;
		}
		for (Hypothesis h : hyps) {
			if (Double.isFinite(h.predIrt) && Double.isFinite(h.retentionTime)) {
				h.predRtError = (float) (Math.abs(h.retentionTime - PredictedSpectrumFeatures.interpolate(bx, by,
						h.predIrt)) / scale);
			}
		}
	}

	static boolean standard(String peptide) {
		for (int i = 0; i < peptide.length(); i++) {
			if ("ACDEFGHIKLMNPQRSTVWY".indexOf(peptide.charAt(i)) < 0) {
				return false;
			}
		}
		return !peptide.isEmpty();
	}
}
