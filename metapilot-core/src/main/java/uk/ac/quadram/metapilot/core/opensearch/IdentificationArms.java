package uk.ac.quadram.metapilot.core.opensearch;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.io.PrintWriter;
import java.nio.charset.StandardCharsets;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.IdentityHashMap;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.Set;
import java.util.zip.GZIPOutputStream;

/**
 * The identification methods compared on one run's searches, so "does the joint identification beat
 * each search alone, their plain union and the cascade" can be answered at matched error rather than
 * at each method's own nominal 1%:
 * <ul>
 * <li>{@code closed_only} - the closed search, scored by a classifier trained on it alone;</li>
 * <li>{@code open_only} - the open search over every spectrum, as the open-only workflow validates it
 * (class {@code open_unexplained} rows are its identification-only tier);</li>
 * <li>{@code cascade} - closed first, open on the spectra the closed stage leaves;</li>
 * <li>{@code joint} - {@link JointIdentification}.</li>
 * </ul>
 * The plain PSM-level union of {@code closed_only} and {@code open_only} is built from these rows by
 * the evaluation script, which also holds the ground truth.
 * <p>
 * Each method's spectrum winners with {@code q <= }{@link #MAX_Q} are written, with the method's own
 * q-value and whether the method accepted the PSM, so the script can trace each method's yield
 * against measured error over a range of thresholds. Rows are written as each method finishes,
 * because the next method re-walks the same PSM objects and overwrites their q-values.
 */
public final class IdentificationArms {

	/** Winners up to this q-value are written, for curves beyond the 1% point. */
	public static final double MAX_Q = 0.05;

	private final PrintWriter rows;
	private final File peptidesFile;
	private final File summaryFile;
	private final Map<String, String> proteins = new HashMap<>();
	private final Map<String, int[]> psmCounts = new LinkedHashMap<>();
	private final Map<String, Set<String>> targetPeptides = new LinkedHashMap<>();
	private final Map<String, Set<String>> decoyPeptides = new LinkedHashMap<>();

	public IdentificationArms(File rowsFile, File peptidesFile, File summaryFile) throws IOException {
		rowsFile.getParentFile().mkdirs();
		this.rows = new PrintWriter(new BufferedWriter(new OutputStreamWriter(
				new GZIPOutputStream(new FileOutputStream(rowsFile), 1 << 16), StandardCharsets.UTF_8), 1 << 20));
		this.rows.println("arm\tfile\tscan\tpeptide\tis_target\tclass\tq_value\taccepted");
		this.peptidesFile = peptidesFile;
		this.summaryFile = summaryFile;
	}

	/** The closed search on its own: its winners and q-values from {@code validateClosed}. */
	public void captureClosedOnly(ClosedOpenCascade.ClosedStage closed, double fdr) {
		for (OpenSearchPsm psm : closed.winners) {
			double q = psm.getQValue();
			if (q <= MAX_Q) {
				row("closed_only", psm, "closed", q, q <= fdr);
			}
		}
	}

	/** The open search on its own, after {@code ModificationValidator.validate} over every spectrum. */
	public void captureOpenOnly(List<OpenSearchPsm> winners, ModificationValidator.Result validation,
			ModificationValidator.Gate gate) {
		Set<OpenSearchPsm> accepted = identitySet(validation.validated);
		Map<Integer, ModificationValidator.MassClass> classes = ModificationValidator.massClasses(gate);
		for (OpenSearchPsm psm : winners) {
			double q = psm.getQValue();
			if (q <= MAX_Q) {
				row("open_only", psm, openClass(psm, gate, classes), q, accepted.contains(psm));
			}
		}
	}

	/** The cascade: closed winners, then open winners of the spectra the closed stage left. */
	public void captureCascade(ClosedOpenCascade.ClosedStage closed, List<OpenSearchPsm> openWinners,
			ModificationValidator.Result cascade, ModificationValidator.Gate gate) {
		Set<OpenSearchPsm> accepted = identitySet(cascade.validated);
		Map<Integer, ModificationValidator.MassClass> classes = ModificationValidator.massClasses(gate);
		for (OpenSearchPsm psm : closed.winners) {
			double q = psm.getQValue();
			if (q <= MAX_Q) {
				row("cascade", psm, "closed", q, accepted.contains(psm));
			}
		}
		for (OpenSearchPsm psm : openWinners) {
			if (closed.acceptedBySpectrum.containsKey(ClosedOpenCascade.spectrumKey(psm))) {
				continue;
			}
			double q = psm.getQValue();
			if (q <= MAX_Q) {
				row("cascade", psm, openClass(psm, gate, classes), q, accepted.contains(psm));
			}
		}
	}

	public void captureJoint(JointIdentification.Result joint) {
		for (int i = 0; i < joint.candidates.size(); i++) {
			if (!joint.winner[i] || !(joint.psmQ[i] <= MAX_Q)) {
				continue;
			}
			JointCandidate c = joint.candidates.get(i);
			row("joint", c.representative(), c.getGroup().name().toLowerCase(Locale.ROOT), joint.psmQ[i],
					joint.kept[i]);
		}
	}

	/** Closes the row table and writes the peptide-protein table and the decoy-estimated summary. */
	public void finish() throws IOException {
		rows.close();
		try (PrintWriter w = new PrintWriter(new BufferedWriter(new OutputStreamWriter(
				new GZIPOutputStream(new FileOutputStream(peptidesFile)), StandardCharsets.UTF_8)))) {
			w.println("peptide\tproteins");
			for (Map.Entry<String, String> e : proteins.entrySet()) {
				w.println(e.getKey() + "\t" + e.getValue());
			}
		}
		Set<String> closedTargets = targetPeptides.getOrDefault("closed_only", Collections.emptySet());
		try (PrintWriter w = new PrintWriter(summaryFile, "UTF-8")) {
			w.println("# Accepted at each method's own 1%: [targets, decoys]. Peptides: base sequence, I = L.");
			w.println("# new_vs_closed_only: accepted target peptides the closed search alone does not accept.");
			w.println("# The ground-truth comparison is the evaluation script's, on identification_arms.tsv.gz.");
			w.println("arm\tpsms_targets\tpsms_decoys\tpeptides_targets\tpeptides_decoys\tpeptide_decoy_ratio\t"
					+ "new_vs_closed_only");
			for (String arm : psmCounts.keySet()) {
				int[] p = psmCounts.get(arm);
				Set<String> t = targetPeptides.getOrDefault(arm, Collections.emptySet());
				Set<String> d = decoyPeptides.getOrDefault(arm, Collections.emptySet());
				int fresh = 0;
				for (String pep : t) {
					if (!closedTargets.contains(pep)) {
						fresh++;
					}
				}
				w.println(arm + "\t" + p[0] + "\t" + p[1] + "\t" + t.size() + "\t" + d.size() + "\t"
						+ String.format(Locale.ROOT, "%.5f", d.size() / (double) Math.max(1, t.size())) + "\t"
						+ (arm.equals("closed_only") ? "" : String.valueOf(fresh)));
			}
		}
	}

	private void row(String arm, OpenSearchPsm psm, String cls, double q, boolean accepted) {
		String peptide = psm.getPeptide();
		rows.println(arm + "\t" + psm.getFileName() + "\t" + psm.getScan() + "\t" + peptide + "\t" + psm.isTarget()
				+ "\t" + cls + "\t" + q + "\t" + accepted);
		if (peptide != null && !proteins.containsKey(peptide)) {
			proteins.put(peptide, joinProteins(psm));
		}
		if (accepted && peptide != null) {
			int[] counts = psmCounts.computeIfAbsent(arm, k -> new int[2]);
			counts[psm.isTarget() ? 0 : 1]++;
			String base = peptide.toUpperCase(Locale.ROOT).replace('I', 'L');
			(psm.isTarget() ? targetPeptides : decoyPeptides).computeIfAbsent(arm, k -> new HashSet<>()).add(base);
		} else {
			psmCounts.computeIfAbsent(arm, k -> new int[2]);
		}
	}

	private static String openClass(OpenSearchPsm psm, ModificationValidator.Gate gate,
			Map<Integer, ModificationValidator.MassClass> classes) {
		int peak = psm.getPeakId();
		if (peak >= 0 && gate.acceptedPeakIds.contains(peak)) {
			ModificationValidator.MassClass mc = classes.get(peak);
			if (mc == ModificationValidator.MassClass.ZERO) {
				return "open_unmodified";
			}
			if (mc == ModificationValidator.MassClass.ESTABLISHED) {
				return "open_accepted_shift";
			}
		}
		return "open_unexplained";
	}

	private static String joinProteins(OpenSearchPsm psm) {
		String[] all = psm.getAllProteins();
		if (all != null && all.length > 0) {
			StringBuilder sb = new StringBuilder();
			for (String p : all) {
				if (p != null && !p.isBlank()) {
					sb.append(sb.length() > 0 ? ";" : "").append(p.trim());
				}
			}
			if (sb.length() > 0) {
				return sb.toString();
			}
		}
		return psm.getProtein() == null ? "" : psm.getProtein();
	}

	private static Set<OpenSearchPsm> identitySet(List<OpenSearchPsm> psms) {
		Set<OpenSearchPsm> set = Collections.newSetFromMap(new IdentityHashMap<>());
		set.addAll(psms);
		return set;
	}
}
