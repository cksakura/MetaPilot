package uk.ac.quadram.metapilot.core.opensearch;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.List;
import java.util.Locale;

/**
 * Whether a delta-mass peak's shift is on the peptide at all, judged from its fragments against what
 * unmodified peptides show at the same mass.
 * <p>
 * Per PSM, {@link ShiftLocaliser} calls a shift fragment-supported when its support reaches
 * {@link ShiftLocaliser#minSupport()}. Some PSMs reach it by chance: dense spectra and a free choice of
 * site give any mass a few matching "shifted" ions. A peak of 26,000 PSMs then yields hundreds of
 * localisations and, with UNIMOD's breadth, explanations, even when nothing at that mass is on the
 * peptide (on 12mix, the +4/+5 Da precursor-error band: {@code Label:13C(3)15N(1)}, kynurenine).
 * <p>
 * The test compares two rates at each peak:
 * <ul>
 * <li>{@code r1}: the share of the peak's accepted winners whose own shift is fragment-supported.
 * Accepted only: below the FDR cut, winners are often wrong peptides, and a wrong peptide in an open
 * search wins by explaining peaks with the free shift, so they inflate {@code r1};</li>
 * <li>{@code r0}: the share of a fixed sample of unmodified zero-peak winners that are
 * fragment-supported when scored at the peak's centroid - the chance rate at that mass - taken at its
 * one-sided 95% upper bound (Wilson), so that a null of 0 to 3 in 2,000 does not read as no chance at
 * all.</li>
 * </ul>
 * {@code r0 / r1} bounds the share of the peak's supported PSMs that are supported by chance. The peak
 * passes when that share is at most {@link #DEFAULT_MAX_CHANCE_SHARE}: at most one localisation in ten
 * expected by chance. A peak that fails is reported as not on the fragments: every winner in it is
 * {@link ShiftLocaliser.State#UNSUPPORTED} and exported as the unmodified peptide. A peak that passes
 * keeps each PSM's own localisation. Neither changes groups, scores or q-values.
 * <p>
 * Measured on the 12mix and faecal entrapment runs (2026-11-07): real modifications 0.003-0.048;
 * -17.027 Da, whose null is high because unmodified spectra carry b-NH3 and y-NH3 ions, 0.052-0.059;
 * the +4/+5 Da band and the +1..+8 ladder 0.079 and above, apart from three peaks of 63-90 PSMs.
 * Nothing here reads a target or decoy label; which winners are accepted comes from the group walks.
 */
public final class PeakFragmentTest {

	/** Largest bounded chance share {@code r0 / r1} at which a peak's shift is on the fragments. */
	public static final double DEFAULT_MAX_CHANCE_SHARE = 0.10;
	/** Unmodified zero-peak winners scored at each peak's mass to estimate the chance rate. */
	public static final int DEFAULT_NULL_SAMPLE = 2000;
	/** One-sided 95% normal quantile for the upper bound on {@code r0}. */
	static final double Z_ONE_SIDED_95 = 1.6448536269514722;

	/** One peak's test. */
	public static final class PeakResult {
		public final int peakId;
		public final double centroid;
		/** Accepted winners in the peak whose localisation was scored. */
		public final int scored;
		/** Of those, fragment-supported. */
		public final int supported;
		/** Unmodified PSMs scored at the centroid. */
		public final int nullScored;
		/** Of those, fragment-supported by chance. */
		public final int nullSupported;
		public final boolean passed;
		public final String reason;

		PeakResult(int peakId, double centroid, int scored, int supported, int nullScored, int nullSupported,
				boolean passed, String reason) {
			this.peakId = peakId;
			this.centroid = centroid;
			this.scored = scored;
			this.supported = supported;
			this.nullScored = nullScored;
			this.nullSupported = nullSupported;
			this.passed = passed;
			this.reason = reason;
		}

		public double observedRate() {
			return scored == 0 ? Double.NaN : supported / (double) scored;
		}

		public double nullRate() {
			return nullScored == 0 ? Double.NaN : nullSupported / (double) nullScored;
		}

		/** One-sided 95% Wilson upper bound on {@link #nullRate()}. */
		public double nullUpperBound() {
			return nullScored == 0 ? Double.NaN : wilsonUpper(nullSupported, nullScored);
		}

		/** {@code r0 upper bound / r1}; infinity when no winner is supported, NaN without a null. */
		public double chanceShare() {
			return PeakFragmentTest.chanceShare(nullUpperBound(), observedRate());
		}
	}

	private PeakFragmentTest() {
	}

	static double chanceShare(double r0, double r1) {
		if (Double.isNaN(r0) || Double.isNaN(r1)) {
			return Double.NaN;
		}
		return r1 <= 0 ? Double.POSITIVE_INFINITY : r0 / r1;
	}

	/** One-sided 95% Wilson score upper bound on {@code k / n}. */
	static double wilsonUpper(int k, int n) {
		double z = Z_ONE_SIDED_95;
		double p = k / (double) n;
		double z2 = z * z;
		return (p + z2 / (2 * n) + z * Math.sqrt(p * (1 - p) / n + z2 / (4.0 * n * n))) / (1 + z2 / n);
	}

	/**
	 * The verdict on one peak. Without a null sample, or without a scored accepted winner, the test
	 * cannot be made and the peak keeps its PSMs' own localisations; the reason says so.
	 */
	static PeakResult decide(int peakId, double centroid, int scored, int supported, int nullScored, int nullSupported,
			double maxChanceShare) {
		if (nullScored == 0) {
			return new PeakResult(peakId, centroid, scored, supported, nullScored, nullSupported, true,
					"not tested: no unmodified PSMs to score");
		}
		if (scored == 0) {
			return new PeakResult(peakId, centroid, scored, supported, nullScored, nullSupported, true,
					"not tested: no scored accepted winners");
		}
		double share = chanceShare(wilsonUpper(nullSupported, nullScored), supported / (double) scored);
		boolean passed = share <= maxChanceShare;
		return new PeakResult(peakId, centroid, scored, supported, nullScored, nullSupported, passed,
				passed ? "on the fragments" : "not on the fragments: chance share above "
						+ String.format(Locale.ROOT, "%.2f", maxChanceShare));
	}

	/** {@code peak_fragment_test.tsv}: one row per tested peak. */
	public static void write(List<PeakResult> results, File out) throws IOException {
		out.getParentFile().mkdirs();
		try (PrintWriter w = new PrintWriter(out, "UTF-8")) {
			w.println("# Per-peak fragment test (PeakFragmentTest): r1 = share of the peak's accepted winners whose shift");
			w.println("# the fragments carry; r0 = the same share for unmodified PSMs scored at the peak's mass (chance),");
			w.println("# r0_upper its one-sided 95% Wilson bound. A peak passes when r0_upper / r1 <= "
					+ DEFAULT_MAX_CHANCE_SHARE + "; failed peaks are exported unmodified.");
			w.println("peak_id\tcentroid_da\twinners_scored\twinners_supported\tr1\tnull_scored\tnull_supported\tr0\t"
					+ "r0_upper\tchance_share\tpassed\treason");
			for (PeakResult r : results) {
				w.println(r.peakId + "\t" + String.format(Locale.ROOT, "%.4f", r.centroid) + "\t" + r.scored + "\t"
						+ r.supported + "\t" + fmt(r.observedRate()) + "\t" + r.nullScored + "\t" + r.nullSupported + "\t"
						+ fmt(r.nullRate()) + "\t" + fmt(r.nullUpperBound()) + "\t" + fmt(r.chanceShare()) + "\t"
						+ r.passed + "\t" + r.reason);
			}
		}
	}

	private static String fmt(double v) {
		return Double.isNaN(v) ? "NaN" : Double.isInfinite(v) ? "Inf" : String.format(Locale.ROOT, "%.4f", v);
	}
}
