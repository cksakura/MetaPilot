package uk.ac.quadram.metapilot.core.opensearch;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.TreeSet;
import java.util.function.Consumer;

import uk.ac.quadram.metapilot.core.variant.MzmlScanReader;

/**
 * Whether a peak is a sequence difference rather than a mass change, from localisation alone.
 * <p>
 * A shift localised to a residue that cannot become another amino acid by that mass is not a
 * substitution. A shift on a residue that can may be either: the fragments are the same, so one
 * PSM cannot be decided. A peak can. For each peak that passed the screen, this takes its rank-1
 * target PSMs that {@link ShiftLocaliser} places on a site, and the share of them on a site where
 * residue plus shift is an amino acid. A peak of wrong peptides sits on those sites; so does a
 * modification whose own residues are substitution sites, deamidation (N/Q, the same mass as N->D
 * and Q->E) above all. The rule cannot tell those two apart and does not try: a withdrawn peak makes
 * no claim, modification or substitution, which is the outcome for anything the fragments cannot
 * decide. The test is arithmetic on the twenty residue masses and the localised site: no list of
 * modifications and no database.
 * <p>
 * What withdrawal does: the peak is no identification group, so its open winners are judged one by
 * one like any unexplained shift. One whose shift the fragments place on a site no substitution
 * fits is identified ({@code OPEN_LOCALISED}); one placed on a substitution site
 * ({@code OPEN_SUBSTITUTION_SITE}) or not placed at all is reported in {@code unexplained_psms.tsv}
 * and not identified. On a substitution site the fragments confirm the mass sits there but not the
 * rest of the peptide: on 12mix (2026-11-08) such PSMs were 11.7% outside the true genomes in the
 * withdrawn peaks (111 PSMs) and 18.8% elsewhere (1,943), against 3.5% for {@code OPEN_LOCALISED}
 * and 2.3% for exact-mass PSMs, and the one-residue variant was a true peptide in almost none of
 * them (5 of 378), so they are wrong peptides rather than variants of the assigned one. Exact-mass
 * PSMs in a withdrawn peak (closed-search isotope-offset matches at 0.984 Da) stay identified as
 * the closed search reported them; only their label is open, and it lists every explanation.
 * <p>
 * Calibration on 12mix against the true genomes, with the localiser before its 13C rule: -14.016 at
 * 0.93 with 20% of PSMs outside them, +14.016 at 0.80 with 20%, +28.031 at 0.80 with 14%, and
 * +40.031 at 0.51 with 1.6%. Below {@link Config#maxSubstitutionSiteShare} of 0.6 the pooled
 * off-truth rate was 2.08%, against 2.06% for the database test and list this replaced, with more
 * PSMs kept; on the four ATCC strains off-strain rates were 10.48, 8.04, 2.54 and 5.44% against
 * 10.36, 7.23, 2.60 and 4.74%. Peak-level, because about half of the +/-14 PSMs have no localised
 * site and are just as wrong.
 * <p>
 * The threshold is not tuned. With the current localiser (round 4, 2026-11-08) the kept peaks of 20
 * or more localised PSMs reach at most 0.50 (12mix +40.031; faecal +1.983 at 0.16) and the withdrawn
 * ones start at 0.72 (faecal +10.021), so any value from 0.51 to 0.72 withdraws the same peaks:
 * 12mix 0.984 (0.80) and +28.031 (0.79); faecal 0.984 (0.87), +/-14.016, +28.031, +30.010, -17.957,
 * -10.021, -23.016, +10.021 and +13.980. Deamidation was at 0.23 under the earlier localiser, which
 * counted 13C peaks as evidence and so placed many 0.984 shifts off N/Q.
 * <p>
 * A peak with fewer than {@link Config#minLocalisedPsms} localised PSMs stays: there is nothing to
 * judge it on, and most such peaks are clean.
 */
public final class PeakSiteShare {

	public static final class Config {
		/** Fewest localised target PSMs a peak needs for its share to count. */
		public int minLocalisedPsms = 20;
		/** A peak at or above this share is withdrawn. */
		public double maxSubstitutionSiteShare = 0.6;
		/** Tolerance for residue plus shift matching another residue, Da. */
		public double massToleranceDa = 0.005;
	}

	/** One peak's counts. */
	public static final class PeakShare {
		public final int peakId;
		public final double centroid;
		public int localised;
		public int onSubstitutionSite;
		public boolean withdrawn;

		PeakShare(int peakId, double centroid) {
			this.peakId = peakId;
			this.centroid = centroid;
		}

		public double share() {
			return localised == 0 ? Double.NaN : onSubstitutionSite / (double) localised;
		}
	}

	public static final class Result {
		public final Map<Integer, PeakShare> peaks = new TreeMap<>();
		/** PSMs of files without readable spectra; these did not count. */
		public int withoutSpectra;

		public Set<Integer> withdrawnPeakIds() {
			Set<Integer> out = new TreeSet<>();
			for (PeakShare p : peaks.values()) {
				if (p.withdrawn) {
					out.add(p.peakId);
				}
			}
			return out;
		}
	}

	private PeakSiteShare() {
	}

	/**
	 * @param psms    open PSMs; only rank-1 targets on {@code peakIds} are read
	 * @param peakIds the non-zero peaks that passed the screen
	 * @param spectra mzML by sample
	 */
	public static Result measure(List<OpenSearchPsm> psms, Set<Integer> peakIds, ShiftLocaliser localiser,
			Map<String, File> spectra, Config cfg, Consumer<String> info) throws IOException {
		Result result = new Result();
		Map<String, List<OpenSearchPsm>> byFile = new LinkedHashMap<>();
		for (OpenSearchPsm psm : psms) {
			if (psm.getHitRank() == 1 && psm.isTarget() && peakIds.contains(psm.getPeakId())) {
				byFile.computeIfAbsent(psm.getFileName(), k -> new ArrayList<>()).add(psm);
				result.peaks.computeIfAbsent(psm.getPeakId(), k -> new PeakShare(k, psm.getPeakCentroid()));
			}
		}
		for (Map.Entry<String, List<OpenSearchPsm>> e : byFile.entrySet()) {
			File mzml = spectra == null ? null : spectra.get(e.getKey());
			if (mzml == null || !mzml.isFile()) {
				result.withoutSpectra += e.getValue().size();
				continue;
			}
			Set<Integer> scans = new HashSet<>();
			for (OpenSearchPsm psm : e.getValue()) {
				scans.add(psm.getScan());
			}
			MzmlScanReader.Result read = MzmlScanReader.read(mzml, scans);
			for (OpenSearchPsm psm : e.getValue()) {
				ShiftLocaliser.Result loc = localiser.localise(psm.getPeptide(), psm.getCharge(), psm.getMassDiff(),
						read.scans.get(psm.getScan()));
				if (loc.state == ShiftLocaliser.State.LOCALISED_EXPLAINED
						|| loc.state == ShiftLocaliser.State.LOCALISED_UNEXPLAINED) {
					count(result.peaks.get(psm.getPeakId()), psm.getPeptide(), loc.site, psm.getMassDiff(), cfg);
				}
			}
		}
		decide(result, cfg);
		if (info != null) {
			info.accept("Site share: " + result.peaks.size() + " peaks tested, " + result.withdrawnPeakIds().size()
					+ " withdrawn as sequence differences (share of localised PSMs on substitution sites >= "
					+ cfg.maxSubstitutionSiteShare + ")"
					+ (result.withoutSpectra > 0 ? "; " + result.withoutSpectra + " PSMs had no readable spectra" : ""));
		}
		return result;
	}

	static void count(PeakShare peak, String peptide, int site, double delta, Config cfg) {
		if (site < 0 || site >= peptide.length()) {
			return;
		}
		peak.localised++;
		if (isSubstitutionSite(peptide.charAt(site), delta, cfg)) {
			peak.onSubstitutionSite++;
		}
	}

	/** Whether {@code residue} plus {@code delta} is another amino acid. */
	static boolean isSubstitutionSite(char residue, double delta, Config cfg) {
		char r = Character.toUpperCase(residue);
		for (SubstitutionMassTable.Substitution s : SubstitutionMassTable.matching(delta, cfg.massToleranceDa)) {
			if (s.from == r || s.from == 'L' && r == 'I' || s.from == 'I' && r == 'L') {
				return true;
			}
		}
		return false;
	}

	static void decide(Result result, Config cfg) {
		for (PeakShare p : result.peaks.values()) {
			p.withdrawn = p.localised >= cfg.minLocalisedPsms && p.share() >= cfg.maxSubstitutionSiteShare;
		}
	}

	public static void write(Result result, File out) throws IOException {
		out.getParentFile().mkdirs();
		try (PrintWriter w = new PrintWriter(out, "UTF-8")) {
			w.println("# Share of a peak's localised rank-1 target PSMs on a site where residue plus shift is an");
			w.println("# amino acid. Peaks that passed the evidence rule only. withdrawn=true: treated as a");
			w.println("# sequence difference, its PSMs make no modification claim.");
			w.println("peak_id\tcentroid_da\tlocalised_psms\ton_substitution_site\tsite_share\twithdrawn");
			for (PeakShare p : result.peaks.values()) {
				w.println(p.peakId + "\t" + p.centroid + "\t" + p.localised + "\t" + p.onSubstitutionSite + "\t"
						+ p.share() + "\t" + p.withdrawn);
			}
		}
	}
}
