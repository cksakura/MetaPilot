package uk.ac.quadram.metapilot.core.unified;

import java.util.Arrays;

import uk.ac.quadram.metapilot.core.opensearch.SubstitutionMassTable;
import uk.ac.quadram.metapilot.core.variant.MzmlScanReader;

/**
 * MetaPilot's own annotation of a spectrum against a hypothesis: every candidate, whichever search
 * proposed it, is scored here by one function, so no search engine score reaches the model.
 * <p>
 * Two ladders per hypothesis:
 * <ul>
 * <li><b>unshifted</b> - the backbone b/y ions (C carbamidomethyl only);</li>
 * <li><b>explained</b> - the ions with the hypothesis's own mass explanation placed: the closed
 * search's variable modifications where it put them, or an open delta on the residue that explains
 * the spectrum best. Equal to the unshifted ladder when there is nothing to place.</li>
 * </ul>
 * Peak picking mirrors the searches' settings: the {@link #TOP_PEAKS} most intense peaks at or
 * above {@link #MIN_RATIO} of the base peak, b/y ions at fragment charge up to
 * {@code min(z - 1, 2)}.
 */
public final class SpectrumAnnotator {

	static final double PROTON = 1.007276466812;
	static final double H2O = 18.010564684;
	static final int TOP_PEAKS = 150;
	static final double MIN_RATIO = 0.01;
	/** Neighbour filtering compares theoretical ions in this m/z range (CONGA). */
	static final double NEIGHBOUR_MIN_MZ = 200;
	static final double NEIGHBOUR_MAX_MZ = 3000;

	private final double tolerancePpm;

	public SpectrumAnnotator(double tolerancePpm) {
		this.tolerancePpm = tolerancePpm;
	}

	/** A spectrum reduced to its picked peaks, sorted by m/z, intensities relative to the base peak. */
	public static final class Spectrum {
		final double[] mz;
		final double[] intensity;
		final double tic;
		final double top50;
		final boolean[] top10;

		Spectrum(double[] mz, double[] intensity, double tic, double top50, boolean[] top10) {
			this.mz = mz;
			this.intensity = intensity;
			this.tic = tic;
			this.top50 = top50;
			this.top10 = top10;
		}

		int size() {
			return mz.length;
		}
	}

	/** Annotation features of one hypothesis. */
	public static final class Annotation {
		float uMatched;
		float uIntTic;
		float uIntTop50;
		float uBCov;
		float uYCov;
		float uBRun;
		float uYRun;
		float uPpmMean;
		float uPpmSd;
		float uTop10;
		/** Hyperscore of the unshifted backbone ladder. */
		float uHyperscore;
		float eMatched;
		float eIntTic;
		float eGain;
		float eSiteMargin;
		float eTiedSites;
		float eTop10;
		float eHyperscore;
		/** 0-based residue carrying an open delta, -1 when none was placed. */
		int bestSite = -1;
		/** Shifted ions gained at the best site over the unshifted ladder (count). */
		int shiftedIonsGained;
	}

	public static Spectrum prepare(MzmlScanReader.Peaks peaks) {
		int n = peaks.mz.length;
		if (n == 0) {
			return new Spectrum(new double[0], new double[0], 0, 0, new boolean[0]);
		}
		Integer[] order = new Integer[n];
		float max = 0f;
		for (int i = 0; i < n; i++) {
			order[i] = i;
			max = Math.max(max, peaks.intensity[i]);
		}
		Arrays.sort(order, (a, b) -> Float.compare(peaks.intensity[b], peaks.intensity[a]));
		int keep = 0;
		while (keep < Math.min(TOP_PEAKS, n) && peaks.intensity[order[keep]] >= MIN_RATIO * max && max > 0) {
			keep++;
		}
		Integer[] kept = Arrays.copyOf(order, keep);
		boolean[] isTop10 = new boolean[n];
		for (int k = 0; k < Math.min(10, keep); k++) {
			isTop10[kept[k]] = true;
		}
		double top50 = 0;
		for (int k = 0; k < Math.min(50, keep); k++) {
			top50 += peaks.intensity[kept[k]] / max;
		}
		Arrays.sort(kept, (a, b) -> Double.compare(peaks.mz[a], peaks.mz[b]));
		double[] mz = new double[keep];
		double[] in = new double[keep];
		boolean[] t10 = new boolean[keep];
		double tic = 0;
		for (int k = 0; k < keep; k++) {
			mz[k] = peaks.mz[kept[k]];
			in[k] = peaks.intensity[kept[k]] / max;
			t10[k] = isTop10[kept[k]];
			tic += in[k];
		}
		return new Spectrum(mz, in, tic, top50, t10);
	}

	/**
	 * Annotates one hypothesis.
	 *
	 * @param mods  closed-search variable modifications (index 0 N-terminus, i residue i), or null
	 * @param delta an open delta to place on the best residue; 0 for none
	 */
	public Annotation annotate(Spectrum s, String peptide, int precursorCharge, double[] mods, double delta) {
		Ladder ladder = new Ladder(peptide, precursorCharge);
		Annotation a = new Annotation();
		int f = ladder.size();
		int[] unshifted = new int[f];
		for (int k = 0; k < f; k++) {
			unshifted[k] = match(s, ladder.mz(k, 0.0));
		}
		fillUnshifted(a, s, ladder, unshifted);

		int[] explained;
		if (mods != null) {
			explained = new int[f];
			for (int k = 0; k < f; k++) {
				explained[k] = match(s, ladder.mzWithMods(k, mods));
			}
		} else if (Math.abs(delta) >= 0.05 && ladder.length > 1) {
			int[] shifted = new int[f];
			for (int k = 0; k < f; k++) {
				shifted[k] = match(s, ladder.mz(k, delta));
			}
			double best = -1;
			double second = -1;
			int bestSite = -1;
			int tied = 0;
			for (int site = 0; site < ladder.length; site++) {
				double score = 0;
				for (int k = 0; k < f; k++) {
					int idx = ladder.contains(k, site) ? shifted[k] : unshifted[k];
					if (idx >= 0) {
						score += s.intensity[idx];
					}
				}
				if (score > best + 1e-12) {
					second = best;
					best = score;
					bestSite = site;
					tied = 1;
				} else if (Math.abs(score - best) <= 1e-12) {
					tied++;
				} else if (score > second) {
					second = score;
				}
			}
			explained = new int[f];
			int gained = 0;
			for (int k = 0; k < f; k++) {
				boolean on = ladder.contains(k, bestSite);
				explained[k] = on ? shifted[k] : unshifted[k];
				if (on && shifted[k] >= 0 && unshifted[k] < 0) {
					gained++;
				}
			}
			a.bestSite = bestSite;
			a.shiftedIonsGained = gained;
			a.eSiteMargin = (float) (s.tic > 0 ? (best - Math.max(0, second)) / s.tic : 0);
			a.eTiedSites = (float) Math.log1p(Math.max(0, tied - 1));
		} else {
			explained = unshifted;
		}
		fillExplained(a, s, ladder, explained);
		return a;
	}

	private void fillUnshifted(Annotation a, Spectrum s, Ladder ladder, int[] idx) {
		boolean[] used = new boolean[s.size()];
		double sum = 0;
		double top10 = 0;
		int matched = 0;
		double ppmSum = 0;
		double ppmSq = 0;
		int ppmN = 0;
		boolean[] b = new boolean[ladder.length];
		boolean[] y = new boolean[ladder.length];
		int nb = 0;
		int ny = 0;
		double ib = 0;
		double iy = 0;
		for (int k = 0; k < idx.length; k++) {
			int i = idx[k];
			if (i < 0) {
				continue;
			}
			(ladder.isY(k) ? y : b)[ladder.ion(k)] = true;
			if (ladder.isY(k)) {
				ny++;
				iy += s.intensity[i];
			} else {
				nb++;
				ib += s.intensity[i];
			}
			double ppm = (s.mz[i] - ladder.mz(k, 0.0)) / ladder.mz(k, 0.0) * 1e6;
			ppmSum += ppm;
			ppmSq += ppm * ppm;
			ppmN++;
			if (!used[i]) {
				used[i] = true;
				matched++;
				sum += s.intensity[i];
				if (s.top10[i]) {
					top10++;
				}
			}
		}
		int ions = Math.max(1, ladder.length - 1);
		a.uMatched = matched;
		a.uIntTic = (float) (s.tic > 0 ? sum / s.tic : 0);
		a.uIntTop50 = (float) (s.top50 > 0 ? Math.min(1.0, sum / s.top50) : 0);
		a.uBCov = (float) count(b) / ions;
		a.uYCov = (float) count(y) / ions;
		a.uBRun = (float) longestRun(b) / ions;
		a.uYRun = (float) longestRun(y) / ions;
		double mean = ppmN == 0 ? 0 : ppmSum / ppmN;
		a.uPpmMean = (float) Math.abs(mean);
		a.uPpmSd = (float) (ppmN < 2 ? tolerancePpm : Math.sqrt(Math.max(0, ppmSq / ppmN - mean * mean)));
		a.uTop10 = (float) (top10 / Math.max(1, Math.min(10, s.size())));
		a.uHyperscore = (float) (lnFactorial(nb) + lnFactorial(ny) + Math.log1p(100 * ib) + Math.log1p(100 * iy));
	}

	private static void fillExplained(Annotation a, Spectrum s, Ladder ladder, int[] idx) {
		boolean[] used = new boolean[s.size()];
		double sum = 0;
		double top10 = 0;
		int matched = 0;
		int nb = 0;
		int ny = 0;
		double ib = 0;
		double iy = 0;
		for (int k = 0; k < idx.length; k++) {
			int i = idx[k];
			if (i < 0) {
				continue;
			}
			if (ladder.isY(k)) {
				ny++;
				iy += s.intensity[i];
			} else {
				nb++;
				ib += s.intensity[i];
			}
			if (!used[i]) {
				used[i] = true;
				matched++;
				sum += s.intensity[i];
				if (s.top10[i]) {
					top10++;
				}
			}
		}
		a.eMatched = matched;
		a.eIntTic = (float) (s.tic > 0 ? sum / s.tic : 0);
		a.eGain = a.eIntTic - a.uIntTic;
		a.eTop10 = (float) (top10 / Math.max(1, Math.min(10, s.size())));
		a.eHyperscore = (float) (lnFactorial(nb) + lnFactorial(ny) + Math.log1p(100 * ib) + Math.log1p(100 * iy));
	}

	/**
	 * Theoretical ions of the explained ladder (charge 1, and 2 when the precursor is 3+) in
	 * {@link #NEIGHBOUR_MIN_MZ}..{@link #NEIGHBOUR_MAX_MZ}, sorted: CONGA's neighbour comparison.
	 */
	static double[] theoreticalIons(String peptide, int precursorCharge, double[] mods, double delta, int site) {
		Ladder ladder = new Ladder(peptide, Math.min(precursorCharge, 3));
		double[] out = new double[ladder.size()];
		int n = 0;
		for (int k = 0; k < ladder.size(); k++) {
			double mz;
			if (mods != null) {
				mz = ladder.mzWithMods(k, mods);
			} else {
				mz = ladder.mz(k, site >= 0 && ladder.contains(k, site) ? delta : 0.0);
			}
			if (mz >= NEIGHBOUR_MIN_MZ && mz <= NEIGHBOUR_MAX_MZ) {
				out[n++] = mz;
			}
		}
		double[] ions = Arrays.copyOf(out, n);
		Arrays.sort(ions);
		return ions;
	}

	/** CONGA's t12 = 2 * shared / (B1 + B2), shared within {@code tolerancePpm}. */
	static double sharedIonFraction(double[] a, double[] b, double tolerancePpm) {
		if (a.length + b.length == 0) {
			return 0;
		}
		int i = 0;
		int j = 0;
		int shared = 0;
		while (i < a.length && j < b.length) {
			double tol = a[i] * tolerancePpm * 1e-6;
			if (Math.abs(a[i] - b[j]) <= tol) {
				shared++;
				i++;
				j++;
			} else if (a[i] < b[j]) {
				i++;
			} else {
				j++;
			}
		}
		return 2.0 * shared / (a.length + b.length);
	}

	/** Index of the most intense peak within tolerance of {@code mz}, or -1. */
	private int match(Spectrum s, double mz) {
		double tol = mz * tolerancePpm * 1e-6;
		int i = Arrays.binarySearch(s.mz, mz - tol);
		if (i < 0) {
			i = -i - 1;
		}
		int best = -1;
		for (; i < s.mz.length && s.mz[i] <= mz + tol; i++) {
			if (best < 0 || s.intensity[i] > s.intensity[best]) {
				best = i;
			}
		}
		return best;
	}

	private static int count(boolean[] v) {
		int c = 0;
		for (boolean x : v) {
			if (x) {
				c++;
			}
		}
		return c;
	}

	private static int longestRun(boolean[] v) {
		int best = 0;
		int run = 0;
		for (int i = 1; i < v.length; i++) {
			run = v[i] ? run + 1 : 0;
			best = Math.max(best, run);
		}
		return best;
	}

	private static final double[] LN_FACT = new double[256];
	static {
		for (int i = 1; i < LN_FACT.length; i++) {
			LN_FACT[i] = LN_FACT[i - 1] + Math.log(i);
		}
	}

	private static double lnFactorial(int n) {
		return n < LN_FACT.length ? LN_FACT[n] : LN_FACT[LN_FACT.length - 1];
	}

	/**
	 * b_i (residues 0..i-1) and y_i (residues n-i..n-1), i = 1..n-1, at fragment charges
	 * 1..max(1, min(z - 1, 2)). Fragment k: ion type, ion number, charge.
	 */
	static final class Ladder {
		final int length;
		final double[] prefix;
		final int maxCharge;

		Ladder(String peptide, int precursorCharge) {
			length = peptide.length();
			prefix = new double[length + 1];
			for (int i = 0; i < length; i++) {
				prefix[i + 1] = prefix[i] + SubstitutionMassTable.residueMass(peptide.charAt(i));
			}
			maxCharge = Math.max(1, Math.min(precursorCharge - 1, 2));
		}

		int size() {
			return 2 * Math.max(0, length - 1) * maxCharge;
		}

		boolean isY(int k) {
			return (k / Math.max(1, length - 1)) % 2 == 1;
		}

		/** Ion number 1..n-1. */
		int ion(int k) {
			return k % Math.max(1, length - 1) + 1;
		}

		int charge(int k) {
			return k / (2 * Math.max(1, length - 1)) + 1;
		}

		/** Whether fragment k holds the residue at 0-based {@code site}. */
		boolean contains(int k, int site) {
			int i = ion(k);
			return isY(k) ? site >= length - i : site <= i - 1;
		}

		double neutral(int k) {
			int i = ion(k);
			return isY(k) ? prefix[length] - prefix[length - i] + H2O : prefix[i];
		}

		double mz(int k, double shift) {
			int z = charge(k);
			return (neutral(k) + shift + z * PROTON) / z;
		}

		double mzWithMods(int k, double[] mods) {
			int i = ion(k);
			double extra = 0;
			if (isY(k)) {
				for (int r = length - i; r < length; r++) {
					extra += r + 1 < mods.length ? mods[r + 1] : 0;
				}
			} else {
				extra += mods.length > 0 ? mods[0] : 0;
				for (int r = 0; r < i; r++) {
					extra += r + 1 < mods.length ? mods[r + 1] : 0;
				}
			}
			return mz(k, extra);
		}
	}
}
