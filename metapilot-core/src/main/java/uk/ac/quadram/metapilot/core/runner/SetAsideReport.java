package uk.ac.quadram.metapilot.core.runner;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.TreeSet;

import uk.ac.quadram.metapilot.core.opensearch.OpenSearchPsm;
import uk.ac.quadram.metapilot.core.opensearch.ShiftLocaliser;

/**
 * What the identification set aside, at the levels the main result is read at: per peptide and per
 * genome. The PSM level is {@code unexplained_psms.tsv}.
 * <p>
 * A set-aside PSM passed its group's 1% walk but is no identification: its shift is a substitution,
 * could be one, or is not accounted for. It is in no peptide, protein, genome or function table, so
 * without these two files a reader of those tables cannot see what was left out, nor that a strain
 * differing from its reference genome loses more PSMs than one that matches it.
 * <ul>
 * <li>{@code set_aside_peptides.tsv}: one row per base peptide, reason and shift mass.</li>
 * <li>{@code set_aside_genomes.tsv}: per genome, set-aside PSMs and peptides by reason beside the
 * genome's identified PSMs and peptides. A peptide shared by several genomes counts for each, on
 * both sides, so shares are comparable between genomes but do not add up to the run's totals.</li>
 * </ul>
 * Targets only. Report-only: nothing reads these files.
 */
final class SetAsideReport {

	static final String PEPTIDES_FILE = "set_aside_peptides.tsv";
	static final String GENOMES_FILE = "set_aside_genomes.tsv";

	/** Why a PSM is set aside, most specific first. */
	enum Reason {
		/** The fragments place the shift on a residue and only an amino acid substitution explains it. */
		SUBSTITUTION,
		/** As above, but a UNIMOD modification fits the same residue: the spectrum cannot decide. */
		MODIFICATION_OR_SUBSTITUTION,
		/** The fragments place the shift on a residue where nothing explains it. */
		SITE_UNEXPLAINED,
		/** The fragments carry the shift but no residue is confident. */
		SITE_AMBIGUOUS,
		/** The fragments do not carry the shift, and UNIMOD has no modification of its mass on this peptide. */
		SHIFT_NOT_ON_FRAGMENTS,
		/** No localisation was made (no readable spectrum). */
		NOT_LOCALISED;

		String label() {
			return name().toLowerCase(Locale.ROOT);
		}
	}

	private SetAsideReport() {
	}

	static Reason reasonOf(OpenSearchPsm psm) {
		ShiftLocaliser.Result loc = psm.getShiftLocalisation();
		if (loc == null || loc.site < 0) {
			return Reason.NOT_LOCALISED;
		}
		switch (loc.state) {
		case LOCALISED_EXPLAINED:
			return loc.explanationClass() == ShiftLocaliser.ExplanationClass.SUBSTITUTION ? Reason.SUBSTITUTION
					: Reason.MODIFICATION_OR_SUBSTITUTION;
		case LOCALISED_UNEXPLAINED:
			return Reason.SITE_UNEXPLAINED;
		case UNSUPPORTED:
			return Reason.SHIFT_NOT_ON_FRAGMENTS;
		default:
			return Reason.SITE_AMBIGUOUS;
		}
	}

	/** Genome of a MAG protein accession ({@code MGYG000002463_01761} -> {@code MGYG000002463}); else the accession. */
	static String genomeOf(String accession) {
		int at = accession.lastIndexOf('_');
		return at > 0 && accession.substring(at + 1).chars().allMatch(Character::isDigit) ? accession.substring(0, at)
				: accession;
	}

	private static Set<String> genomesOf(OpenSearchPsm psm, String decoyTag) {
		Set<String> out = new TreeSet<>();
		String[] all = psm.getAllProteins();
		List<String> accessions = new ArrayList<>();
		if (all != null) {
			Collections.addAll(accessions, all);
		}
		if (accessions.isEmpty() && psm.getProtein() != null) {
			accessions.add(psm.getProtein());
		}
		for (String a : accessions) {
			if (a != null && !a.isBlank() && (decoyTag == null || !a.trim().startsWith(decoyTag))) {
				out.add(genomeOf(a.trim()));
			}
		}
		return out;
	}

	private static final class PeptideRow {
		String peptide;
		Reason reason;
		final List<Double> masses = new ArrayList<>();
		final Set<String> explanations = new TreeSet<>();
		final Set<Integer> sites = new TreeSet<>();
		final Set<String> samples = new TreeSet<>();
		final Set<String> genomes = new TreeSet<>();
		final Set<String> groups = new TreeSet<>();
		int psms;
		double bestQ = Double.POSITIVE_INFINITY;
	}

	private static final class GenomeRow {
		final int[] psms = new int[Reason.values().length];
		final List<Set<String>> peptides = new ArrayList<>();
		int identifiedPsms;
		final Set<String> identifiedPeptides = new TreeSet<>();

		GenomeRow() {
			for (int i = 0; i < Reason.values().length; i++) {
				peptides.add(new TreeSet<>());
			}
		}

		int setAsidePsms() {
			int n = 0;
			for (int v : psms) {
				n += v;
			}
			return n;
		}

		int setAsidePeptides() {
			Set<String> all = new TreeSet<>();
			for (Set<String> s : peptides) {
				all.addAll(s);
			}
			return all.size();
		}
	}

	/**
	 * @param setAside   the PSMs of {@code unexplained_psms.tsv}
	 * @param identified the PSMs of {@code validated_psms.tsv}
	 * @return rows written to the peptide table
	 */
	static int write(List<OpenSearchPsm> setAside, List<OpenSearchPsm> identified, File dir, String decoyTag)
			throws IOException {
		Map<String, PeptideRow> peptides = new LinkedHashMap<>();
		Map<String, GenomeRow> genomes = new TreeMap<>();
		if (setAside != null) {
			for (OpenSearchPsm psm : setAside) {
				if (psm == null || !psm.isTarget() || psm.getPeptide() == null) {
					continue;
				}
				Reason reason = reasonOf(psm);
				double mass = psm.getPeakId() >= 0 && Double.isFinite(psm.getNormalizedMassDiff())
						? psm.getNormalizedMassDiff() : psm.getMassDiff();
				String key = psm.getPeptide() + '\t' + reason + '\t' + String.format(Locale.ROOT, "%.2f", mass);
				PeptideRow row = peptides.computeIfAbsent(key, k -> new PeptideRow());
				row.peptide = psm.getPeptide();
				row.reason = reason;
				row.masses.add(mass);
				row.psms++;
				row.samples.add(String.valueOf(psm.getFileName()));
				row.bestQ = Math.min(row.bestQ, psm.getQValue());
				if (psm.getIdentificationGroup() != null) {
					row.groups.add(psm.getIdentificationGroup());
				}
				ShiftLocaliser.Result loc = psm.getShiftLocalisation();
				if (loc != null && psm.getLocalisedShiftSite() >= 0) {
					row.sites.add(psm.getLocalisedShiftSite() + 1);
					row.explanations.addAll(loc.explanations);
				}
				Set<String> psmGenomes = genomesOf(psm, decoyTag);
				row.genomes.addAll(psmGenomes);
				for (String g : psmGenomes) {
					GenomeRow gr = genomes.computeIfAbsent(g, k -> new GenomeRow());
					gr.psms[reason.ordinal()]++;
					gr.peptides.get(reason.ordinal()).add(psm.getPeptide());
				}
			}
		}
		if (identified != null) {
			for (OpenSearchPsm psm : identified) {
				if (psm == null || !psm.isTarget() || psm.getPeptide() == null) {
					continue;
				}
				for (String g : genomesOf(psm, decoyTag)) {
					GenomeRow gr = genomes.computeIfAbsent(g, k -> new GenomeRow());
					gr.identifiedPsms++;
					gr.identifiedPeptides.add(psm.getPeptide());
				}
			}
		}
		dir.mkdirs();
		try (PrintWriter w = new PrintWriter(new File(dir, PEPTIDES_FILE), "UTF-8")) {
			w.println("# PSMs that passed their group's FDR walk but are not identifications, per base peptide, reason");
			w.println("# and shift mass. They are in no peptide, protein, genome or function table. Targets only.");
			w.println("# reason: substitution = only an amino acid substitution explains the shift at its site;");
			w.println("#   modification_or_substitution = a UNIMOD modification fits the same residue;");
			w.println("#   site_unexplained = the fragments place the shift where nothing explains it;");
			w.println("#   site_ambiguous = the fragments carry the shift but no residue is confident;");
			w.println("#   shift_not_on_fragments = not on the fragments, and no UNIMOD modification of that mass fits");
			w.println("#   the peptide; not_localised = no readable spectrum.");
			w.println("peptide\treason\tshift_mass_da\tshift_sites\tshift_explanation\tid_groups\tpsms\tsamples_observed\t"
					+ "best_q_value\tgenomes\tsamples");
			for (PeptideRow r : peptides.values()) {
				List<Double> m = new ArrayList<>(r.masses);
				Collections.sort(m);
				List<String> sites = new ArrayList<>();
				for (Integer s : r.sites) {
					sites.add(String.valueOf(s));
				}
				w.println(r.peptide + "\t" + r.reason.label() + "\t"
						+ String.format(Locale.ROOT, "%.4f", m.get(m.size() / 2)) + "\t" + String.join(";", sites) + "\t"
						+ String.join("|", r.explanations) + "\t" + String.join(";", r.groups) + "\t" + r.psms + "\t"
						+ r.samples.size() + "\t" + String.format(Locale.ROOT, "%.6g", r.bestQ) + "\t"
						+ String.join(";", r.genomes) + "\t" + String.join(";", r.samples));
			}
		}
		try (PrintWriter w = new PrintWriter(new File(dir, GENOMES_FILE), "UTF-8")) {
			w.println("# Per genome: PSMs and base peptides set aside (set_aside_peptides.tsv), by reason, beside the");
			w.println("# genome's identified PSMs and peptides (validated_psms.tsv). A peptide shared by several");
			w.println("# genomes counts for each of them on both sides. set_aside_psm_share = set aside / (set aside +");
			w.println("# identified): a genome whose strain in the sample differs from the reference loses more.");
			StringBuilder header = new StringBuilder("genome\tidentified_psms\tidentified_peptides\tset_aside_psms\t"
					+ "set_aside_peptides\tset_aside_psm_share");
			for (Reason reason : Reason.values()) {
				header.append('\t').append(reason.label()).append("_psms\t").append(reason.label()).append("_peptides");
			}
			w.println(header);
			for (Map.Entry<String, GenomeRow> e : genomes.entrySet()) {
				GenomeRow g = e.getValue();
				int aside = g.setAsidePsms();
				StringBuilder row = new StringBuilder(e.getKey()).append('\t').append(g.identifiedPsms).append('\t')
						.append(g.identifiedPeptides.size()).append('\t').append(aside).append('\t')
						.append(g.setAsidePeptides()).append('\t')
						.append(aside + g.identifiedPsms == 0 ? ""
								: String.format(Locale.ROOT, "%.4f", aside / (double) (aside + g.identifiedPsms)));
				for (Reason reason : Reason.values()) {
					row.append('\t').append(g.psms[reason.ordinal()]).append('\t')
							.append(g.peptides.get(reason.ordinal()).size());
				}
				w.println(row);
			}
		}
		return peptides.size();
	}
}
