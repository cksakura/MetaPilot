package uk.ac.quadram.metapilot.core.unified;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.io.PrintWriter;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.function.Consumer;
import java.util.stream.IntStream;
import java.util.zip.GZIPOutputStream;

import uk.ac.quadram.metapilot.core.opensearch.DeltaMassPeakFitter;
import uk.ac.quadram.metapilot.core.opensearch.OpenSearchPsm;
import uk.ac.quadram.metapilot.core.variant.MzmlScanReader;

/**
 * Unified identification, phase 1 of {@code docs/developer/unified-rescoring-batchwise-fdr-design.md}.
 * <ol>
 * <li>Delta-mass peaks detected label-blind, filtered (F2-F6), mapped afterwards
 * ({@link DeltaMassBatches}).</li>
 * <li>Per spectrum, one candidate set from the closed and the open search; every candidate
 * annotated by MetaPilot ({@link SpectrumAnnotator}), optionally scored against DIA-NN's backbone
 * prediction ({@link DiannBackbone}); neighbours filtered (CONGA, tau 5%).</li>
 * <li>A listwise model, one query per spectrum with a "none" option ({@link ListwiseGbdt}), trained
 * semi-supervised with the RESET decoy split and folds by spectrum.</li>
 * <li>Each spectrum's winner among its reportable candidates, scored by its log-odds against every
 * alternative including "none"; Group-walk over the batches at PSM and at peptide level, with the
 * RESET estimate.</li>
 * </ol>
 * No search-engine score is a feature: the searches only propose candidates.
 */
public final class UnifiedIdentification {

	static final double ISOTOPE = 1.0033548;

	public static final class Config {
		public double fragmentTolerancePpm = 20.0;
		public double fdr = 0.01;
		public double trainQ = 0.01;
		public int folds = 3;
		public int iterations = 3;
		public int groupWalkK = 40;
		public boolean reportUnexplained = false;
		/**
		 * Shared-ion fraction above which two candidates are neighbours. CONGA's 5% made unrelated
		 * short peptides neighbours on one coincidental ion (a 7-mer has ~12 ions in range): on 12mix,
		 * 36,777 of 43,422 closed-accepted spectra that lost their exact candidate this way had an
		 * unrelated replacement (2026-10-20). 50% keeps the filter to real hitchhikers (I/L, one or two
		 * substitutions, nested sequences) and leaves the rest to the model.
		 */
		public double neighbourTau = 0.5;
		public int maxTrainingQueries = 150_000;
		/** Diagnostic: candidates from the closed search only. */
		public boolean closedOnly = false;
		/**
		 * Diagnostic: the closed classifier's feature set ({@code ClosedOpenCascade.CLOSED_FEATURES}:
		 * MSFragger scores and the same-peptide/protein counts) instead of MetaPilot's, to check the
		 * training and FDR machinery against the reference. Not part of the method.
		 */
		public boolean engineFeatures = false;
		/**
		 * RESET (train on half the decoys, estimate with the other half). Off: Percolator-style, every
		 * decoy trains and every decoy estimates.
		 */
		public boolean reset = true;
		/** Ablation: feature names left out of the model. */
		public Set<String> dropFeatures = new HashSet<>();
		/**
		 * Stage-wise, exact-first identification ({@code docs/developer/stagewise-exact-first-identification-design.md}):
		 * exact-mass candidates of both searches decided first, shifted candidates only on the spectra
		 * stage 1 left. Off: one competition over every candidate (the unified arm).
		 */
		public boolean stagewise = false;
		/**
		 * Recurrence counts (same peptide, same protein) as PSM-score features. Off by default: with them
		 * the RT-consistency check doubled its RT-inconsistent share at the same decoy estimate on 12mix
		 * (3.2% -> 6.2%). Recurrence stays in the peptide-level model. An ablation arm until validated by
		 * entrapment.
		 */
		public boolean psmRecurrence = false;
		public ListwiseGbdt.Params gbdt = new ListwiseGbdt.Params();
		public DeltaMassBatches.Config batches = new DeltaMassBatches.Config();
		public DeltaMassPeakFitter.Config fitter = new DeltaMassPeakFitter.Config();
		/** DIA-NN 1.8.x executable; null turns the prediction block off. */
		public File diann;
		public File workDir;
		public int threads = Runtime.getRuntime().availableProcessors();
	}

	/** Per spectrum winner. */
	static final class Winner {
		final Hypothesis hyp;
		final double score;
		final double utility;
		Winner(Hypothesis hyp, double score, double utility) {
			this.hyp = hyp;
			this.score = score;
			this.utility = utility;
		}
	}

	private final Config cfg;
	private final Consumer<String> log;
	private final List<String> specFile = new ArrayList<>();
	private final List<Integer> specScan = new ArrayList<>();
	private String[] closedTop;
	private String[] openTop;
	private final List<Hypothesis> hyps = new ArrayList<>();
	private DeltaMassBatches batches;
	private final List<String> summary = new ArrayList<>();
	/** Peptide-unit table written for diagnosis; set from the output folder. */
	private File unitTable;

	private UnifiedIdentification(Config cfg, Consumer<String> log) {
		this.cfg = cfg;
		this.log = log;
	}

	public static List<String> run(List<OpenSearchPsm> open, List<OpenSearchPsm> closed, Map<String, File> mzml,
			Config cfg, File outDir, Consumer<String> log) throws IOException {
		UnifiedIdentification u = new UnifiedIdentification(cfg, log);
		outDir.mkdirs();
		u.unitTable = new File(outDir, "unified_peptide_units.tsv.gz");
		long t0 = System.currentTimeMillis();
		u.buildCandidates(open, closed);
		u.detectPeaks(open);
		u.annotate(mzml);
		u.filterPeaksAndExplain(open);
		u.time("candidates, peaks, annotation", t0);
		long t1 = System.currentTimeMillis();
		u.scoreAndWalk(outDir);
		u.time("model and FDR", t1);
		return u.summary;
	}

	private void time(String what, long t0) {
		String line = String.format(Locale.ROOT, "%s: %.1f min", what, (System.currentTimeMillis() - t0) / 60000.0);
		log.accept("Unified: " + line);
		summary.add("time_" + what.replace(", ", "_").replace(' ', '_') + "\t"
				+ String.format(Locale.ROOT, "%.1f", (System.currentTimeMillis() - t0) / 60000.0));
	}

	// ------------------------------------------------------------ 1. candidate sets

	private void buildCandidates(List<OpenSearchPsm> open, List<OpenSearchPsm> closed) {
		Map<String, Integer> ids = new HashMap<>();
		Map<String, Hypothesis> byKey = new HashMap<>(Math.max(16, open.size() + closed.size()));
		Map<Integer, String> closedTopMap = new HashMap<>();
		Map<Integer, String> openTopMap = new HashMap<>();
		int conflicts = 0;
		for (int pass = cfg.closedOnly ? 1 : 0; pass < 2; pass++) {
			List<OpenSearchPsm> rows = pass == 0 ? open : closed;
			for (int i = 0; i < rows.size(); i++) {
				OpenSearchPsm psm = rows.get(i);
				if (psm.getPeptide() == null || psm.getFileName() == null) {
					continue;
				}
				String sk = psm.getFileName() + '\t' + psm.getScan();
				Integer sid = ids.get(sk);
				if (sid == null) {
					sid = ids.size();
					ids.put(sk, sid);
					specFile.add(psm.getFileName());
					specScan.add(psm.getScan());
				}
				// A peptidoform: closed rows by sequence and variable modifications, open rows by sequence;
				// an open row joins a closed one only when they explain the same mass (mergeEquivalent).
				String key = sid + "\t" + psm.getPeptide() + (pass == 0 ? "\tO" : "\tC" + modSignature(psm.getClosedModifications()));
				Hypothesis h = byKey.get(key);
				if (h == null) {
					h = new Hypothesis(sid, psm.getPeptide(), psm.isTarget(), psm.getCharge(),
							psm.getPrecursorNeutralMass(), psm.getRetentionTime(), psm.getNumTolTerm(),
							psm.getNumMissedCleavages(), psm.getAllProteins() != null && psm.getAllProteins().length > 0
									? psm.getAllProteins() : new String[] { psm.getProtein() });
					byKey.put(key, h);
					hyps.add(h);
				} else if (h.target != psm.isTarget()) {
					conflicts++;
				}
				if (pass == 0) {
					h.foundByOpen = true;
					if (psm.getHitRank() < h.openRank) {
						h.openRank = psm.getHitRank();
						h.openMassDiff = psm.getMassDiff();
						h.openFitIndex = i;
					}
					if (psm.getHitRank() == 1) {
						openTopMap.put(sid, Hypothesis.unit(psm.getPeptide()));
					}
				} else {
					h.foundByClosed = true;
					if (psm.getHitRank() < h.closedRank) {
						h.closedRank = psm.getHitRank();
						h.closedMassDiff = psm.getMassDiff();
						h.closedMods = psm.getClosedModifications();
						h.protein = psm.getProtein();
						if (cfg.engineFeatures) {
							h.engine = new float[] { (float) psm.getHyperscore(), (float) psm.getNextscore(),
									(float) (psm.getHyperscore() - psm.getNextscore()),
									(float) -Math.log10(Math.max(psm.getExpect(), 1e-300)), psm.getNumMatchedIons(),
									psm.getTotNumIons(), (float) psm.getMatchedIonFraction(),
									psm.getClosedModificationCount() };
						}
					}
					if (psm.getHitRank() == 1) {
						closedTopMap.put(sid, Hypothesis.unit(psm.getPeptide()));
					}
				}
			}
		}
		closedTop = new String[ids.size()];
		openTop = new String[ids.size()];
		closedTopMap.forEach((k, v) -> closedTop[k] = v);
		openTopMap.forEach((k, v) -> openTop[k] = v);
		int merged = mergeEquivalent(hyps);
		for (Hypothesis h : hyps) {
			// RESET roles per decoy peptide: every PSM of a decoy peptide, in every stage, shares the role,
			// so the peptide level keeps the training/estimating separation. Without RESET every decoy
			// trains (and every decoy also estimates, see walkAndWrite).
			h.resetTraining = !h.target && (!cfg.reset || (stableHash("decoy:" + Hypothesis.unit(h.peptide)) & 1L) == 0L);
		}
		summary.add("open_candidates_merged_into_equivalent_closed\t" + merged);
		if (cfg.engineFeatures) {
			// The closed classifier's aggregate counts, over the closed rank-1 rows as it computes them.
			Map<String, Integer> pep = new HashMap<>();
			Map<String, Integer> prot = new HashMap<>();
			for (Hypothesis h : hyps) {
				if (h.closedRank == 1) {
					pep.merge(h.peptide, 1, Integer::sum);
					if (h.protein != null) {
						prot.merge(h.protein, 1, Integer::sum);
					}
				}
			}
			for (Hypothesis h : hyps) {
				if (h.engine != null) {
					float[] e = Arrays.copyOf(h.engine, h.engine.length + 2);
					e[h.engine.length] = (float) Math.log1p(pep.getOrDefault(h.peptide, 0));
					e[h.engine.length + 1] = (float) Math.log1p(h.protein == null ? 0 : prot.getOrDefault(h.protein, 0));
					h.engine = e;
				}
			}
		}
		log.accept(String.format(Locale.ROOT,
				"Unified: %,d open + %,d closed rows -> %,d spectra, %,d candidates (%,d target/decoy label conflicts)",
				open.size(), closed.size(), ids.size(), hyps.size(), conflicts));
		summary.add("spectra\t" + ids.size());
		summary.add("candidates\t" + hyps.size());
	}

	/** Variable modifications as a key: position and mass to 0.01 Da; empty when unmodified. */
	static String modSignature(double[] mods) {
		if (mods == null) {
			return "";
		}
		StringBuilder sb = new StringBuilder();
		for (int i = 0; i < mods.length; i++) {
			if (Math.abs(mods[i]) > 1e-3) {
				sb.append(i).append(':').append(Math.round(mods[i] * 100)).append(';');
			}
		}
		return sb.toString();
	}

	/**
	 * Merges an open-only hypothesis into a closed hypothesis of the same spectrum and sequence when
	 * both explain the same mass: open delta = the closed candidate's modification mass + its
	 * precursor offset (isotope), within 0.02 Da. The merged hypothesis is "found by both". Otherwise
	 * the open candidate stays a peptidoform of its own (the same backbone with another explanation).
	 */
	static int mergeEquivalent(List<Hypothesis> hyps) {
		Map<String, List<Hypothesis>> closedBy = new HashMap<>();
		for (Hypothesis h : hyps) {
			if (h.foundByClosed) {
				closedBy.computeIfAbsent(h.spectrum + "\t" + h.peptide, k -> new ArrayList<>()).add(h);
			}
		}
		int merged = 0;
		List<Hypothesis> kept = new ArrayList<>(hyps.size());
		for (Hypothesis h : hyps) {
			if (h.foundByClosed || !h.foundByOpen) {
				kept.add(h);
				continue;
			}
			Hypothesis into = null;
			for (Hypothesis c : closedBy.getOrDefault(h.spectrum + "\t" + h.peptide, List.of())) {
				double closedMass = modMass(c.closedMods) + (Double.isFinite(c.closedMassDiff) ? c.closedMassDiff : 0);
				if (Double.isFinite(h.openMassDiff) && Math.abs(h.openMassDiff - closedMass) <= 0.02) {
					into = c;
					break;
				}
			}
			if (into == null) {
				kept.add(h);
			} else {
				into.foundByOpen = true;
				if (h.openRank < into.openRank) {
					into.openRank = h.openRank;
					into.openMassDiff = h.openMassDiff;
					into.openFitIndex = h.openFitIndex;
				}
				merged++;
			}
		}
		hyps.clear();
		hyps.addAll(kept);
		return merged;
	}

	static double modMass(double[] mods) {
		double m = 0;
		if (mods != null) {
			for (double v : mods) {
				m += v;
			}
		}
		return m;
	}

	// ------------------------------------------------------------ 2. peaks

	private void detectPeaks(List<OpenSearchPsm> open) {
		int n = open.size();
		double[] delta = new double[n];
		boolean[] rank1 = new boolean[n];
		double[] weight = new double[n];
		String[] files = new String[n];
		Map<String, Integer> perSpectrum = new HashMap<>();
		for (int i = 0; i < n; i++) {
			OpenSearchPsm p = open.get(i);
			delta[i] = p.getMassDiff();
			rank1[i] = p.getHitRank() == 1;
			files[i] = p.getFileName();
			if (rank1[i]) {
				perSpectrum.merge(p.getFileName() + '\t' + p.getScan(), 1, Integer::sum);
			}
		}
		for (int i = 0; i < n; i++) {
			if (rank1[i]) {
				OpenSearchPsm p = open.get(i);
				weight[i] = 1.0 / perSpectrum.get(p.getFileName() + '\t' + p.getScan());
			}
		}
		batches = DeltaMassBatches.detect(delta, rank1, weight, files, cfg.fitter, log);
		for (Hypothesis h : hyps) {
			if (h.foundByClosed) {
				h.explanation = Hypothesis.Explanation.EXACT;
				int k = Double.isFinite(h.closedMassDiff) ? (int) Math.round(h.closedMassDiff / ISOTOPE) : 0;
				h.isotopeOffset = Math.max(0, Math.min(2, k));
				h.fragmentDelta = 0;
			} else {
				int i = h.openFitIndex;
				h.fragmentDelta = batches.fit.normalizedMassDiff[i];
				h.isotopeOffset = batches.fit.isotopeOffset[i];
				if (batches.isZero(i)) {
					h.explanation = Hypothesis.Explanation.EXACT;
					h.fragmentDelta = 0;
				}
			}
		}
	}

	// ------------------------------------------------------------ 3. annotation

	private void annotate(Map<String, File> mzml) throws IOException {
		SpectrumAnnotator annotator = new SpectrumAnnotator(cfg.fragmentTolerancePpm);
		Map<String, List<Hypothesis>> byFile = new TreeMap<>();
		for (Hypothesis h : hyps) {
			byFile.computeIfAbsent(specFile.get(h.spectrum), k -> new ArrayList<>()).add(h);
		}
		int missingSpectra = 0;
		int neighbours = 0;
		for (Map.Entry<String, List<Hypothesis>> e : byFile.entrySet()) {
			File file = mzml.get(e.getKey());
			List<Hypothesis> fileHyps = e.getValue();
			if (file == null) {
				log.accept("Unified: WARNING no mzML for " + e.getKey() + "; its " + fileHyps.size()
						+ " candidates are left out (targets and decoys alike)");
				for (Hypothesis h : fileHyps) {
					h.neighbour = true;
				}
				missingSpectra += fileHyps.size();
				continue;
			}
			Set<Integer> scans = new HashSet<>();
			for (Hypothesis h : fileHyps) {
				scans.add(specScan.get(h.spectrum));
			}
			long t = System.currentTimeMillis();
			MzmlScanReader.Result read = MzmlScanReader.read(file, scans);
			Map<Integer, SpectrumAnnotator.Spectrum> spectra = new HashMap<>(scans.size() * 2);
			for (Hypothesis h : fileHyps) {
				if (!spectra.containsKey(h.spectrum)) {
					MzmlScanReader.Peaks p = read.scans.get(specScan.get(h.spectrum));
					if (p != null) {
						spectra.put(h.spectrum, SpectrumAnnotator.prepare(p));
					}
				}
			}
			read = null;
			fileHyps.parallelStream().forEach(h -> {
				SpectrumAnnotator.Spectrum s = spectra.get(h.spectrum);
				if (s != null && h.peptide.length() > 1) {
					h.annotation = annotator.annotate(s, h.peptide, h.charge, h.foundByClosed ? h.closedMods : null,
							h.foundByClosed ? 0.0 : h.fragmentDelta);
				}
			});
			for (Hypothesis h : fileHyps) {
				if (h.annotation == null) {
					h.neighbour = true;
					missingSpectra++;
				}
			}
			neighbours += filterNeighbours(fileHyps);
			if (cfg.diann != null && cfg.diann.isFile()) {
				DiannBackbone.scoreRun(fileHyps, spectra, cfg.diann, new File(cfg.workDir, "diann_" + e.getKey()),
						cfg.threads, cfg.fragmentTolerancePpm, log, m -> log.accept("Unified: WARNING " + m));
			}
			log.accept(String.format(Locale.ROOT, "Unified: annotated %,d candidates of %s (%,d spectra) in %.1f s",
					fileHyps.size(), e.getKey(), spectra.size(), (System.currentTimeMillis() - t) / 1000.0));
		}
		summary.add("candidates_without_spectrum\t" + missingSpectra);
		summary.add("neighbours_removed\t" + neighbours);
		summary.add("prediction_block\t" + (cfg.diann != null && cfg.diann.isFile() ? cfg.diann : "off (no DIA-NN)"));
		int predicted = 0;
		int withRt = 0;
		int scored = 0;
		for (Hypothesis h : hyps) {
			if (!h.neighbour) {
				scored++;
				if (!Float.isNaN(h.predEntropy)) {
					predicted++;
				}
				if (!Float.isNaN(h.predRtError)) {
					withRt++;
				}
			}
		}
		summary.add("candidates_with_prediction_[ms2,rt,of]\t[" + predicted + ", " + withRt + ", " + scored + "]");
		if (cfg.diann != null && cfg.diann.isFile() && predicted == 0) {
			log.accept("Unified: WARNING DIA-NN was configured but produced no predictions; see the DIA-NN messages above");
		}
	}

	/**
	 * The label-free score that ranks a spectrum's candidates against each other (neighbour filter,
	 * context features, first ranking): an exact candidate's own ladder (its variable modifications
	 * placed), an open candidate's <em>unshifted</em> backbone ladder. The open shift earns no fragment
	 * credit here, as in CONGA's XCorr comparison of narrow and open hits. Crediting it let an open
	 * variant of the exact peptide score a little higher by placing a shift on its fragments: on 12mix
	 * 53,832 closed-accepted spectra lost their exact candidate to such a neighbour (2026-10-20).
	 */
	static float rankScore(Hypothesis h) {
		return h.foundByClosed ? h.annotation.eHyperscore : h.annotation.uHyperscore;
	}

	/**
	 * CONGA's neighbour filter within each spectrum, ranked by {@link #rankScore}. In stage-wise mode
	 * within each stage only: an exact and a shifted hypothesis never remove each other, because they
	 * are decided in different stages.
	 */
	private int filterNeighbours(List<Hypothesis> fileHyps) {
		Map<Long, List<Hypothesis>> bySpectrum = new HashMap<>();
		for (Hypothesis h : fileHyps) {
			if (!h.neighbour) {
				long stage = cfg.stagewise && h.explanation == Hypothesis.Explanation.EXACT ? 1L : 0L;
				bySpectrum.computeIfAbsent(h.spectrum * 2L + stage, k -> new ArrayList<>()).add(h);
			}
		}
		int[] removed = new int[1];
		bySpectrum.values().parallelStream().forEach(list -> {
			list.sort((a, b) -> Float.compare(rankScore(b), rankScore(a)));
			List<double[]> keptIons = new ArrayList<>();
			for (Hypothesis h : list) {
				double[] ions = SpectrumAnnotator.theoreticalIons(h.peptide, h.charge,
						h.foundByClosed ? h.closedMods : null, h.foundByClosed ? 0.0 : h.fragmentDelta,
						h.annotation.bestSite);
				boolean neighbour = false;
				for (double[] k : keptIons) {
					if (SpectrumAnnotator.sharedIonFraction(ions, k, cfg.fragmentTolerancePpm) > cfg.neighbourTau) {
						neighbour = true;
						break;
					}
				}
				if (neighbour) {
					h.neighbour = true;
					synchronized (removed) {
						removed[0]++;
					}
				} else {
					keptIons.add(ions);
				}
			}
		});
		return removed[0];
	}

	// ------------------------------------------------------------ 4. peak filters, explanations

	private void filterPeaksAndExplain(List<OpenSearchPsm> open) {
		int n = open.size();
		boolean[] rank1 = new boolean[n];
		int[] gained = new int[n];
		Arrays.fill(gained, Integer.MIN_VALUE);
		for (int i = 0; i < n; i++) {
			rank1[i] = open.get(i).getHitRank() == 1;
		}
		for (Hypothesis h : hyps) {
			if (!h.foundByClosed && h.openRank == 1 && h.annotation != null) {
				gained[h.openFitIndex] = Math.abs(h.fragmentDelta) >= 0.05 ? h.annotation.shiftedIonsGained : 0;
			}
		}
		batches.filter(rank1, gained, batches.fit.normalizedMassDiff, cfg.batches, log);
		int[] counts = new int[3];
		for (Hypothesis h : hyps) {
			if (h.explanation == null) {
				int fam = batches.familyOf(h.openFitIndex);
				int b = fam >= 0 ? batches.batchOf(fam) : -1;
				if (b >= 0) {
					h.explanation = Hypothesis.Explanation.PEAK;
					h.peakId = b;
				} else {
					h.explanation = Hypothesis.Explanation.UNEXPLAINED;
				}
			}
			if (!h.neighbour) {
				counts[h.explanation.ordinal()]++;
			}
		}
		summary.add("candidates_exact\t" + counts[0]);
		summary.add("candidates_peak\t" + counts[1]);
		summary.add("candidates_unexplained\t" + counts[2]);
	}

	// ------------------------------------------------------------ 5-7. features, model, FDR

	static final String[] FEATURES = { "charge", "peptide_length", "missed_cleavages", "tol_term", "log_precursor_mass",
			"u_matched", "u_int_tic", "u_int_top50", "u_b_cov", "u_y_cov", "u_b_run", "u_y_run", "u_ppm_mean",
			"u_ppm_sd", "u_top10", "u_hyperscore", "e_matched", "e_int_tic", "e_gain", "e_site_margin", "e_tied_sites",
			"e_top10",
			"e_hyperscore", "expl_error_sigma", "isotope_offset", "log_abs_delta", "is_exact", "is_peak",
			"is_unexplained", "log_batch_size", "found_by_closed", "found_by_open", "found_by_both",
			"peak_log_enrichment", "peak_support_rate", "peak_log_sigma", "log_local_delta_density", "ctx_rank",
			"ctx_margin_all", "ctx_margin_exact", "ctx_margin_nonexact", "ctx_n_close", "ctx_search_agree",
			"ctx_rival_exact", "ctx_rival_peak", "ctx_rival_unexplained", "ctx_rival_found_by_both",
			"ctx_rival_peak_log_enrichment", "ctx_rival_peak_support_rate", "ctx_n_candidates", "log_pep_psms_any",
			"log_pep_psms_top", "pep_run_fraction_top", "log_prot_psms_top", "log_prot_peptides_top",
			"ctx_shifted_margin", "ctx_shifted_rival_log_enrichment", "ctx_shifted_rival_support",
			"ctx_shifted_rival_found_by_both", "ctx_has_shifted_rival" };
	/** Recurrence features: in the PSM score only with {@link Config#psmRecurrence}. */
	static final String[] RECURRENCE_FEATURES = { "log_pep_psms_any", "log_pep_psms_top", "pep_run_fraction_top",
			"log_prot_psms_top", "log_prot_peptides_top" };
	static final String[] PRED_FEATURES = { "p_entropy", "p_pearson", "p_observed", "p_rt_error", "p_missing" };
	/** Diagnostic only: the closed classifier's inputs ({@code ClosedOpenCascade.CLOSED_FEATURES}). */
	static final String[] ENGINE_FEATURES = { "charge", "log_precursor_mass", "peptide_length", "tol_term",
			"missed_cleavages", "hyperscore", "nextscore", "hyperscore_minus_nextscore", "neg_log10_expect",
			"num_matched_ions", "tot_num_ions", "matched_ion_fraction", "closed_variable_modifications",
			"log_same_peptide_count", "log_same_protein_count", "normalized_mass_diff_ppm", "isotope_offset" };

	/** Diagnostic: rows of {@link #ENGINE_FEATURES}; candidates without a closed row get zeros. */
	private float[][] engineMatrix(List<Hypothesis> rows) {
		float[][] x = new float[rows.size()][];
		for (int r = 0; r < rows.size(); r++) {
			Hypothesis h = rows.get(r);
			float[] f = new float[ENGINE_FEATURES.length];
			f[0] = h.charge;
			f[1] = (float) Math.log(Math.max(1, h.precursorMass));
			f[2] = h.peptide.length();
			f[3] = Math.max(0, h.numTolTerm);
			f[4] = Math.max(0, h.missedCleavages);
			if (h.engine != null) {
				System.arraycopy(h.engine, 0, f, 5, h.engine.length);
				f[15] = (float) (Double.isFinite(h.closedMassDiff) ? closedPpm(h) : 0);
			}
			f[16] = h.isotopeOffset;
			x[r] = f;
		}
		return x;
	}

	/** Mann-Whitney AUC of target over decoy scores; NaN without both. */
	static double auc(List<double[]> scored) {
		scored.sort((a, b) -> Double.compare(a[0], b[0]));
		double rankSum = 0;
		long nt = 0;
		long nd = 0;
		int i = 0;
		while (i < scored.size()) {
			int j = i;
			while (j + 1 < scored.size() && scored.get(j + 1)[0] == scored.get(i)[0]) {
				j++;
			}
			double mid = (i + j) / 2.0 + 1;
			for (int k = i; k <= j; k++) {
				if (scored.get(k)[1] == 1) {
					rankSum += mid;
					nt++;
				} else {
					nd++;
				}
			}
			i = j + 1;
		}
		return nt == 0 || nd == 0 ? Double.NaN : (rankSum - nt * (nt + 1) / 2.0) / ((double) nt * nd);
	}

	/** Target-vs-decoy AUC of the winners' scores over the spectra {@code fold[s] == f} (or != f). */
	private static double winnerAuc(List<Hypothesis> rows, int[] qStart, double[] util, int[] fold, int f,
			boolean heldOut) {
		List<double[]> scored = new ArrayList<>();
		for (int s = 0; s < qStart.length - 1; s++) {
			if ((fold[s] == f) != heldOut || qStart[s] == qStart[s + 1]) {
				continue;
			}
			Winner w = winner(rows, qStart, util, s);
			scored.add(new double[] { w.score, w.hyp.target ? 1 : 0 });
		}
		return auc(scored);
	}

	private void scoreAndWalk(File outDir) throws IOException {
		boolean pred = false;
		for (Hypothesis h : hyps) {
			if (!Float.isNaN(h.predEntropy)) {
				pred = true;
				break;
			}
		}
		String[] names = cfg.engineFeatures ? ENGINE_FEATURES : pred ? concat(FEATURES, PRED_FEATURES) : FEATURES;
		int[] monotone = cfg.engineFeatures ? new int[names.length] : monotone(names);

		// Rows grouped by spectrum.
		List<List<Hypothesis>> bySpectrum = new ArrayList<>();
		for (int s = 0; s < specFile.size(); s++) {
			bySpectrum.add(new ArrayList<>());
		}
		for (Hypothesis h : hyps) {
			if (!h.neighbour) {
				bySpectrum.get(h.spectrum).add(h);
			}
		}
		int spectra = bySpectrum.size();
		int[] qStart = new int[spectra + 1];
		List<Hypothesis> rows = new ArrayList<>();
		for (int s = 0; s < spectra; s++) {
			qStart[s] = rows.size();
			rows.addAll(bySpectrum.get(s));
		}
		qStart[spectra] = rows.size();
		bySpectrum = null;
		int n = rows.size();
		float[][] x = cfg.engineFeatures ? engineMatrix(rows) : features(rows, qStart, pred, names.length);
		if (!cfg.engineFeatures && !cfg.psmRecurrence) {
			cfg.dropFeatures.addAll(Arrays.asList(RECURRENCE_FEATURES));
		}
		summary.add("psm_recurrence\t" + (cfg.engineFeatures ? "engine feature set" : cfg.psmRecurrence));
		summary.add("mode_stagewise\t" + cfg.stagewise);
		if (!cfg.dropFeatures.isEmpty()) {
			// Ablation: remove the named columns.
			List<Integer> keep = new ArrayList<>();
			for (int k = 0; k < names.length; k++) {
				if (!cfg.dropFeatures.contains(names[k])) {
					keep.add(k);
				}
			}
			String[] kn = new String[keep.size()];
			int[] km = new int[keep.size()];
			for (int k = 0; k < keep.size(); k++) {
				kn[k] = names[keep.get(k)];
				km[k] = monotone[keep.get(k)];
			}
			for (int r = 0; r < x.length; r++) {
				float[] f = new float[keep.size()];
				for (int k = 0; k < keep.size(); k++) {
					f[k] = x[r][keep.get(k)];
				}
				x[r] = f;
			}
			summary.add("features_dropped\t" + String.join(",", cfg.dropFeatures));
			names = kn;
			monotone = km;
		}
		log.accept(String.format(Locale.ROOT, "Unified: %,d candidate rows x %d features (%s)", n, names.length,
				pred ? "with DIA-NN block" : "no prediction block"));
		summary.add("features\t" + String.join(",", names));

		int[] fold = new int[spectra];
		for (int s = 0; s < spectra; s++) {
			fold[s] = spectrumFold(s);
		}
		if (cfg.stagewise) {
			stagewise(rows, qStart, x, names, monotone, outDir);
			return;
		}
		double[] util = fitQueries(rows, qStart, x, fold, monotone, names, "");
		// Final winners and scores.
		Winner[] winners = new Winner[spectra];
		IntStream.range(0, spectra).parallel().forEach(s -> winners[s] = winner(rows, qStart, util, s));
		walkAndWrite(winners, outDir);
	}

	private int spectrumFold(int s) {
		return (int) Math.floorMod(stableHash(specFile.get(s) + ":" + specScan.get(s)), (long) cfg.folds);
	}

	/**
	 * Semi-supervised, cross-validated listwise training over the queries {@code qStart} (one per
	 * spectrum), folds {@code fold} per query. Returns every row's utility from the model of the fold
	 * that held its query out. Diagnostics go to the summary with {@code label} as prefix.
	 */
	private double[] fitQueries(List<Hypothesis> rows, int[] qStart, float[][] x, int[] fold, int[] monotone,
			String[] names, String label) {
		int spectra = qStart.length - 1;
		int n = rows.size();
		double[] util = new double[n];
		double[] base = new double[n];
		for (int r = 0; r < n; r++) {
			base[r] = rankScore(rows.get(r));
		}
		int[] splitCounts = new int[names.length];
		for (int f = 0; f < cfg.folds; f++) {
			double[] working = base.clone();
			ListwiseGbdt.Model model = null;
			for (int it = 0; it < cfg.iterations; it++) {
				TrainingSet ts = trainingSet(rows, qStart, x, working, fold, f);
				log.accept(String.format(Locale.ROOT,
						"Unified: %sfold %d/%d iteration %d - %,d positive spectra, %,d none spectra (%,d rows; %,d "
								+ "positive and %,d none spectra before the cap)",
						label, f + 1, cfg.folds, it + 1, ts.positives, ts.negatives, ts.x.length, ts.positivesBeforeCap,
						ts.negativesBeforeCap));
				summary.add(String.format(Locale.ROOT,
						"%sfold_%d_iteration_%d_[positive,none]_spectra_[used,before_cap]\t[%d, %d, %d, %d]", label,
						f + 1, it + 1, ts.positives, ts.negatives, ts.positivesBeforeCap, ts.negativesBeforeCap));
				if (ts.positives < 50 || ts.negatives < 50) {
					log.accept("Unified: WARNING too few training spectra; keeping the last model");
					break;
				}
				model = ListwiseGbdt.train(ts.x, ts.qStart, ts.positive, monotone, cfg.gbdt);
				final ListwiseGbdt.Model m = model;
				final int trainingOf = f;
				IntStream.range(0, spectra).parallel().filter(s -> fold[s] != trainingOf).forEach(s -> {
					for (int r = qStart[s]; r < qStart[s + 1]; r++) {
						working[r] = m.utility(x[r]);
					}
				});
			}
			final ListwiseGbdt.Model m = model;
			final int heldOut = f;
			IntStream.range(0, spectra).parallel().filter(s -> fold[s] == heldOut).forEach(s -> {
				for (int r = qStart[s]; r < qStart[s + 1]; r++) {
					util[r] = m == null ? base[r] : m.utility(x[r]);
				}
			});
			double trainAuc = winnerAuc(rows, qStart, working, fold, f, false);
			double heldAuc = winnerAuc(rows, qStart, util, fold, f, true);
			summary.add(String.format(Locale.ROOT, "%sfold_%d_winner_auc_[training,held_out]\t[%.4f, %.4f]", label,
					f + 1, trainAuc, heldAuc));
			if (model == null) {
				log.accept("Unified: WARNING fold " + (f + 1) + " has no model; held-out spectra keep the annotation score");
			} else {
				for (ListwiseGbdt.Tree t : model.trees) {
					for (int feat : t.feature) {
						if (feat >= 0) {
							splitCounts[feat]++;
						}
					}
				}
			}
		}
		StringBuilder sc = new StringBuilder();
		Integer[] order = new Integer[names.length];
		for (int i = 0; i < order.length; i++) {
			order[i] = i;
		}
		Arrays.sort(order, (a, b) -> Integer.compare(splitCounts[b], splitCounts[a]));
		for (int i : order) {
			sc.append(sc.length() > 0 ? "," : "").append(names[i]).append(':').append(splitCounts[i]);
		}
		summary.add(label + "tree_splits_per_feature\t" + sc);
		return util;
	}

	private static final class TrainingSet {
		float[][] x;
		int[] qStart;
		int[] positive;
		int positives;
		int negatives;
		int positivesBeforeCap;
		int negativesBeforeCap;
	}

	/**
	 * Training queries of the partition {@code fold != heldOut}: spectra whose winner is a pseudo-target
	 * (target or estimating decoy, RESET) at q <= trainQ within its explanation type become positive
	 * queries; spectra whose winner is a training decoy become "none" queries.
	 */
	private TrainingSet trainingSet(List<Hypothesis> rows, int[] qStart, float[][] x, double[] util, int[] fold,
			int heldOut) {
		int spectra = qStart.length - 1;
		List<Integer> positiveSpectra = new ArrayList<>();
		List<Integer> negativeSpectra = new ArrayList<>();
		Map<Hypothesis.Explanation, List<double[]>> byType = new HashMap<>();
		Winner[] w = new Winner[spectra];
		for (int s = 0; s < spectra; s++) {
			if (fold[s] == heldOut || qStart[s] == qStart[s + 1]) {
				continue;
			}
			w[s] = winner(rows, qStart, util, s);
			Hypothesis h = w[s].hyp;
			if (!h.target && h.resetTraining) {
				negativeSpectra.add(s);
			}
			// {score, spectrum, isTrainingDecoy}
			byType.computeIfAbsent(h.explanation, k -> new ArrayList<>())
					.add(new double[] { w[s].score, s, !h.target && h.resetTraining ? 1 : 0 });
		}
		for (List<double[]> list : byType.values()) {
			list.sort((a, b) -> Double.compare(b[0], a[0]));
			int m = list.size();
			double[] fdr = new double[m];
			int pseudo = 0;
			int training = 0;
			for (int i = 0; i < m; i++) {
				if (list.get(i)[2] == 1) {
					training++;
				} else {
					pseudo++;
				}
				// RESET: about 3 false pseudo-targets per training decoy (false targets 2, estimating
				// decoys 1). Without RESET the pseudo-targets are the targets and every decoy counts once.
				fdr[i] = pseudo == 0 ? 1.0 : Math.min(1.0, (cfg.reset ? 3.0 : 1.0) * training / pseudo);
			}
			double best = 1.0;
			for (int i = m - 1; i >= 0; i--) {
				best = Math.min(best, fdr[i]);
				fdr[i] = best;
			}
			for (int i = 0; i < m; i++) {
				if (list.get(i)[2] == 0 && fdr[i] <= cfg.trainQ) {
					positiveSpectra.add((int) list.get(i)[1]);
				}
			}
		}
		// Cap by whole spectra, ratio preserved, deterministic.
		List<Integer> chosen = new ArrayList<>();
		List<Boolean> isPos = new ArrayList<>();
		int total = positiveSpectra.size() + negativeSpectra.size();
		double keep = total > cfg.maxTrainingQueries ? cfg.maxTrainingQueries / (double) total : 1.0;
		for (int s : positiveSpectra) {
			if (keep >= 1.0 || unit(s) < keep) {
				chosen.add(s);
				isPos.add(true);
			}
		}
		for (int s : negativeSpectra) {
			if (keep >= 1.0 || unit(s) < keep) {
				chosen.add(s);
				isPos.add(false);
			}
		}
		TrainingSet ts = new TrainingSet();
		ts.positivesBeforeCap = positiveSpectra.size();
		ts.negativesBeforeCap = negativeSpectra.size();
		int rowCount = 0;
		for (int s : chosen) {
			rowCount += qStart[s + 1] - qStart[s];
		}
		ts.x = new float[rowCount][];
		ts.qStart = new int[chosen.size() + 1];
		ts.positive = new int[chosen.size()];
		int r = 0;
		for (int q = 0; q < chosen.size(); q++) {
			int s = chosen.get(q);
			ts.qStart[q] = r;
			ts.positive[q] = -1;
			for (int row = qStart[s]; row < qStart[s + 1]; row++) {
				if (isPos.get(q) && rows.get(row) == w[s].hyp) {
					ts.positive[q] = r;
				}
				ts.x[r++] = x[row];
			}
			if (isPos.get(q)) {
				ts.positives++;
			} else {
				ts.negatives++;
			}
		}
		ts.qStart[chosen.size()] = r;
		return ts;
	}

	private static double unit(int s) {
		return (stableHash("q" + s) >>> 11) / (double) (1L << 53);
	}

	/**
	 * The spectrum's winner: highest utility among its reportable candidates (exact, retained peak),
	 * an unexplained one only when there is no reportable candidate. Its score is its log-odds
	 * against every alternative, any explanation or "none".
	 */
	static Winner winner(List<Hypothesis> rows, int[] qStart, double[] util, int s) {
		int from = qStart[s];
		int to = qStart[s + 1];
		if (from == to) {
			return null;
		}
		int best = -1;
		for (int r = from; r < to; r++) {
			if (rows.get(r).explanation != Hypothesis.Explanation.UNEXPLAINED
					&& (best < 0 || util[r] > util[best])) {
				best = r;
			}
		}
		if (best < 0) {
			for (int r = from; r < to; r++) {
				if (best < 0 || util[r] > util[best]) {
					best = r;
				}
			}
		}
		// Alternatives: "none" and the reportable candidates. An unexplained rival's random delta is
		// evidence against itself, not against an explained winner, so it is not subtracted here; the
		// model still sees it through the rival features. A winner that is itself unexplained (no
		// reportable candidate) is weighed against the other unexplained ones.
		boolean winnerReportable = rows.get(best).explanation != Hypothesis.Explanation.UNEXPLAINED;
		double max = 0;
		for (int r = from; r < to; r++) {
			if (r != best && (!winnerReportable || rows.get(r).explanation != Hypothesis.Explanation.UNEXPLAINED)) {
				max = Math.max(max, util[r]);
			}
		}
		double z = Math.exp(-max);
		for (int r = from; r < to; r++) {
			if (r != best && (!winnerReportable || rows.get(r).explanation != Hypothesis.Explanation.UNEXPLAINED)) {
				z += Math.exp(util[r] - max);
			}
		}
		double score = util[best] - (max + Math.log(z));
		return new Winner(rows.get(best), score, util[best]);
	}

	/** Recurrence counts, computed once from the candidate sets ({@link #countRecurrence}). */
	private final Map<String, Integer> pepAny = new HashMap<>();
	private final Map<String, Integer> pepTop = new HashMap<>();
	private final Map<String, Integer> pepRuns = new HashMap<>();
	private final Map<String, Integer> protTopPsms = new HashMap<>();
	private final Map<String, Integer> protTopPeptides = new HashMap<>();
	private int runCount;

	/**
	 * Recurrence of each peptide and protein across the dataset, from the candidate sets alone: how
	 * many spectra propose a peptide (any candidate, either search), how many have it as their top
	 * candidate by the label-free {@link #rankScore}, in how many runs, and for each protein its
	 * top-candidate PSMs and distinct top-candidate peptides. Targets and decoys ({@code rev_}) are
	 * counted the same way; no label and no classifier score is read, so nothing feeds back across
	 * the training iterations.
	 */
	private void countRecurrence(List<Hypothesis> rows, int[] qStart) {
		Map<String, Set<String>> runsOf = new HashMap<>();
		Map<String, Set<String>> peptidesOf = new HashMap<>();
		Set<String> runs = new HashSet<>();
		for (int s = 0; s < qStart.length - 1; s++) {
			int from = qStart[s];
			int to = qStart[s + 1];
			if (from == to) {
				continue;
			}
			Hypothesis top = rows.get(from);
			Set<String> seen = new HashSet<>();
			for (int r = from; r < to; r++) {
				Hypothesis h = rows.get(r);
				if (seen.add(Hypothesis.unit(h.peptide))) {
					pepAny.merge(Hypothesis.unit(h.peptide), 1, Integer::sum);
				}
				if (rankScore(h) > rankScore(top)) {
					top = h;
				}
			}
			String unit = Hypothesis.unit(top.peptide);
			String run = specFile.get(top.spectrum);
			runs.add(run);
			pepTop.merge(unit, 1, Integer::sum);
			runsOf.computeIfAbsent(unit, k -> new HashSet<>()).add(run);
			for (String p : top.proteins) {
				if (p != null && !p.isBlank()) {
					protTopPsms.merge(p, 1, Integer::sum);
					peptidesOf.computeIfAbsent(p, k -> new HashSet<>()).add(unit);
				}
			}
		}
		runsOf.forEach((k, v) -> pepRuns.put(k, v.size()));
		peptidesOf.forEach((k, v) -> protTopPeptides.put(k, v.size()));
		runCount = runs.size();
		summary.add("recurrence_counts\t" + pepTop.size() + " top-candidate peptides, " + protTopPsms.size()
				+ " proteins, " + runCount + " runs");
	}

	private float[][] features(List<Hypothesis> rows, int[] qStart, boolean pred, int width) {
		countRecurrence(rows, qStart);
		int n = rows.size();
		float[][] x = new float[n][];
		// Precursor error scale per run (closed rank-1, label-blind).
		Map<String, double[]> ppmStats = closedPpmStats();
		Map<String, Integer> batchSize = new HashMap<>();
		for (Hypothesis h : rows) {
			batchSize.merge(batchKey(h), 1, Integer::sum);
		}
		List<DeltaMassBatches.Candidate> retained = new ArrayList<>();
		for (DeltaMassBatches.Candidate c : batches.candidates.values()) {
			if (c.batch == c.id) {
				retained.add(c);
			}
		}
		DeltaMassBatches.Candidate zero = batches.fit.zeroPeakId >= 0 ? batches.candidate(batches.fit.zeroPeakId) : null;
		int spectra = qStart.length - 1;
		IntStream.range(0, spectra).parallel().forEach(s -> {
			int from = qStart[s];
			int to = qStart[s + 1];
			double best = Double.NEGATIVE_INFINITY;
			for (int r = from; r < to; r++) {
				best = Math.max(best, rankScore(rows.get(r)));
			}
			boolean agree = closedTop[s] != null && closedTop[s].equals(openTop[s]);
			for (int r = from; r < to; r++) {
				Hypothesis h = rows.get(r);
				SpectrumAnnotator.Annotation a = h.annotation;
				double own = rankScore(h);
				int rank = 1;
				double otherAll = Double.NEGATIVE_INFINITY;
				double otherExact = Double.NEGATIVE_INFINITY;
				double otherNon = Double.NEGATIVE_INFINITY;
				Hypothesis rival = null;
				Hypothesis shiftedRival = null;
				int close = 0;
				for (int k = from; k < to; k++) {
					if (k == r) {
						continue;
					}
					Hypothesis o = rows.get(k);
					double v = rankScore(o);
					if (v > own) {
						rank++;
					}
					if (v >= 0.9 * best) {
						close++;
					}
					if (v > otherAll) {
						otherAll = v;
						rival = o;
					}
					if (o.explanation == Hypothesis.Explanation.EXACT) {
						otherExact = Math.max(otherExact, v);
					} else {
						if (v > otherNon) {
							shiftedRival = o;
						}
						otherNon = Math.max(otherNon, v);
					}
				}
				float[] f = new float[width];
				int i = 0;
				f[i++] = h.charge;
				f[i++] = h.peptide.length();
				f[i++] = Math.max(0, h.missedCleavages);
				f[i++] = Math.max(0, h.numTolTerm);
				f[i++] = (float) Math.log(Math.max(1, h.precursorMass));
				f[i++] = a.uMatched;
				f[i++] = a.uIntTic;
				f[i++] = a.uIntTop50;
				f[i++] = a.uBCov;
				f[i++] = a.uYCov;
				f[i++] = a.uBRun;
				f[i++] = a.uYRun;
				f[i++] = a.uPpmMean;
				f[i++] = a.uPpmSd;
				f[i++] = a.uTop10;
				f[i++] = a.uHyperscore;
				f[i++] = a.eMatched;
				f[i++] = a.eIntTic;
				f[i++] = a.eGain;
				f[i++] = a.eSiteMargin;
				f[i++] = a.eTiedSites;
				f[i++] = a.eTop10;
				f[i++] = a.eHyperscore;
				f[i++] = (float) explanationError(h, ppmStats, zero, retained);
				f[i++] = h.isotopeOffset;
				f[i++] = (float) Math.log1p(Math.abs(h.foundByClosed ? 0 : h.fragmentDelta));
				f[i++] = h.explanation == Hypothesis.Explanation.EXACT ? 1 : 0;
				f[i++] = h.explanation == Hypothesis.Explanation.PEAK ? 1 : 0;
				f[i++] = h.explanation == Hypothesis.Explanation.UNEXPLAINED ? 1 : 0;
				f[i++] = (float) Math.log1p(batchSize.getOrDefault(batchKey(h), 0));
				f[i++] = h.foundByClosed ? 1 : 0;
				f[i++] = h.foundByOpen ? 1 : 0;
				f[i++] = h.foundByClosed && h.foundByOpen ? 1 : 0;
				double[] peak = peakQuality(h, zero);
				f[i++] = (float) peak[0];
				f[i++] = (float) peak[1];
				f[i++] = (float) peak[2];
				f[i++] = (float) peak[3];
				f[i++] = rank;
				f[i++] = (float) (own - (otherAll == Double.NEGATIVE_INFINITY ? 0 : otherAll));
				f[i++] = (float) (own - (otherExact == Double.NEGATIVE_INFINITY ? 0 : otherExact));
				f[i++] = (float) (own - (otherNon == Double.NEGATIVE_INFINITY ? 0 : otherNon));
				f[i++] = close;
				f[i++] = agree ? 1 : 0;
				f[i++] = rival != null && rival.explanation == Hypothesis.Explanation.EXACT ? 1 : 0;
				f[i++] = rival != null && rival.explanation == Hypothesis.Explanation.PEAK ? 1 : 0;
				f[i++] = rival != null && rival.explanation == Hypothesis.Explanation.UNEXPLAINED ? 1 : 0;
				f[i++] = rival != null && rival.foundByClosed && rival.foundByOpen ? 1 : 0;
				double[] rivalPeak = rival == null ? new double[] { 0, 0, 0, 0 } : peakQuality(rival, zero);
				f[i++] = (float) rivalPeak[0];
				f[i++] = (float) rivalPeak[1];
				f[i++] = to - from;
				String unit = Hypothesis.unit(h.peptide);
				f[i++] = (float) Math.log1p(pepAny.getOrDefault(unit, 0));
				f[i++] = (float) Math.log1p(pepTop.getOrDefault(unit, 0));
				f[i++] = pepRuns.getOrDefault(unit, 0) / (float) Math.max(1, runCount);
				int protPsms = 0;
				int protPeps = 0;
				for (String p : h.proteins) {
					if (p != null) {
						protPsms = Math.max(protPsms, protTopPsms.getOrDefault(p, 0));
						protPeps = Math.max(protPeps, protTopPeptides.getOrDefault(p, 0));
					}
				}
				f[i++] = (float) Math.log1p(protPsms);
				f[i++] = (float) Math.log1p(protPeps);
				// What the open search offers against this candidate as a shifted explanation.
				double[] sq = shiftedRival == null ? new double[] { 0, 0, 0, 0 } : peakQuality(shiftedRival, zero);
				f[i++] = (float) (own - (shiftedRival == null ? 0 : rankScore(shiftedRival)));
				f[i++] = (float) sq[0];
				f[i++] = (float) sq[1];
				f[i++] = shiftedRival != null && shiftedRival.foundByClosed && shiftedRival.foundByOpen ? 1 : 0;
				f[i++] = shiftedRival == null ? 0 : 1;
				if (pred) {
					boolean missing = Float.isNaN(h.predEntropy);
					f[i++] = missing ? 0 : h.predEntropy;
					f[i++] = missing ? 0 : h.predPearson;
					f[i++] = missing ? 0 : h.predObserved;
					f[i++] = Float.isNaN(h.predRtError) ? 0 : Math.min(50f, h.predRtError);
					f[i++] = missing ? 1 : 0;
				}
				x[r] = f;
			}
		});
		return x;
	}

	/**
	 * How trustworthy a candidate's mass explanation is, from the peak filters' label-blind
	 * measurements: {log enrichment over the flanks, fragment-support rate, log sigma (Da), log local
	 * delta density}. Exact candidates carry the zero peak's values with support 1 (no shift to
	 * support); a peak candidate its own family root's (not the peak it was merged into); an
	 * unexplained one enrichment 0, the off-peak background support rate and its local density.
	 */
	private double[] peakQuality(Hypothesis h, DeltaMassBatches.Candidate zero) {
		double density = h.foundByClosed || h.openFitIndex < 0 ? 0
				: Math.log1p(Math.max(0, batches.fit.localDensity[h.openFitIndex]));
		if (h.explanation == Hypothesis.Explanation.EXACT) {
			double e = zero == null || Double.isNaN(zero.enrichment) ? 0 : Math.log1p(zero.enrichment);
			return new double[] { e, 1.0, zero == null ? 0 : Math.log(zero.sigma), density };
		}
		if (h.explanation == Hypothesis.Explanation.PEAK) {
			DeltaMassBatches.Candidate own = batches.candidate(batches.familyOf(h.openFitIndex));
			DeltaMassBatches.Candidate c = own != null ? own : batches.candidate(h.peakId);
			double support = c.supportTested == 0 ? batches.backgroundSupportRate
					: c.supportSuccesses / (double) c.supportTested;
			return new double[] { Math.log1p(Double.isNaN(c.enrichment) ? 0 : c.enrichment), support,
					Math.log(c.sigma), density };
		}
		return new double[] { 0, Double.isNaN(batches.backgroundSupportRate) ? 0 : batches.backgroundSupportRate,
				0, density };
	}

	private String batchKey(Hypothesis h) {
		return h.explanation == Hypothesis.Explanation.PEAK ? "PEAK_" + h.peakId : h.explanation.name();
	}

	/** Per run: median and robust SD of the closed rank-1 precursor error after the isotope offset (ppm). */
	private Map<String, double[]> closedPpmStats() {
		Map<String, List<Double>> byFile = new HashMap<>();
		for (Hypothesis h : hyps) {
			if (h.foundByClosed && h.closedRank == 1 && Double.isFinite(h.closedMassDiff)) {
				byFile.computeIfAbsent(specFile.get(h.spectrum), k -> new ArrayList<>()).add(closedPpm(h));
			}
		}
		Map<String, double[]> out = new HashMap<>();
		for (Map.Entry<String, List<Double>> e : byFile.entrySet()) {
			double[] v = e.getValue().stream().mapToDouble(Double::doubleValue).toArray();
			double med = median(v);
			double[] dev = new double[v.length];
			for (int i = 0; i < v.length; i++) {
				dev[i] = Math.abs(v[i] - med);
			}
			out.put(e.getKey(), new double[] { med, Math.max(0.5, 1.4826 * median(dev)) });
		}
		return out;
	}

	private static double closedPpm(Hypothesis h) {
		return (h.closedMassDiff - h.isotopeOffset * ISOTOPE) / Math.max(1, h.precursorMass) * 1e6;
	}

	private double explanationError(Hypothesis h, Map<String, double[]> ppmStats, DeltaMassBatches.Candidate zero,
			List<DeltaMassBatches.Candidate> retained) {
		double cap = 50;
		if (h.foundByClosed) {
			double[] st = ppmStats.get(specFile.get(h.spectrum));
			if (st == null || !Double.isFinite(h.closedMassDiff)) {
				return 0;
			}
			return Math.min(cap, Math.abs(closedPpm(h) - st[0]) / st[1]);
		}
		double d = batches.fit.normalizedMassDiff[h.openFitIndex];
		if (h.explanation == Hypothesis.Explanation.EXACT) {
			return zero == null ? 0 : Math.min(cap, Math.abs(d - zero.centroid) / zero.sigma);
		}
		if (h.explanation == Hypothesis.Explanation.PEAK) {
			DeltaMassBatches.Candidate own = batches.candidate(batches.familyOf(h.openFitIndex));
			DeltaMassBatches.Candidate c = own != null ? own : batches.candidate(h.peakId);
			return Math.min(cap, Math.abs(d - c.centroid) / c.sigma);
		}
		double best = cap;
		for (DeltaMassBatches.Candidate c : retained) {
			best = Math.min(best, Math.abs(d - c.centroid) / c.sigma);
		}
		return best;
	}

	// ------------------------------------------------------------ stage-wise, exact first

	/** One stage's queries: a subset of each spectrum's candidates, with features and folds. */
	private static final class StageQueries {
		final List<Hypothesis> rows = new ArrayList<>();
		final List<float[]> x = new ArrayList<>();
		final List<Integer> qStart = new ArrayList<>();
		final List<Integer> spectrum = new ArrayList<>();
		final List<Integer> fold = new ArrayList<>();

		int[] qStartArray() {
			int[] q = new int[qStart.size() + 1];
			for (int i = 0; i < qStart.size(); i++) {
				q[i] = qStart.get(i);
			}
			q[qStart.size()] = rows.size();
			return q;
		}

		int[] foldArray() {
			return fold.stream().mapToInt(Integer::intValue).toArray();
		}

		float[][] xArray() {
			return x.toArray(new float[0][]);
		}
	}

	private StageQueries stageQueries(List<Hypothesis> rows, int[] qStart, float[][] x,
			java.util.function.IntPredicate spectrumOk, java.util.function.Predicate<Hypothesis> candidateOk,
			java.util.function.IntFunction<float[]> extra) {
		StageQueries out = new StageQueries();
		for (int s = 0; s < qStart.length - 1; s++) {
			if (!spectrumOk.test(s)) {
				continue;
			}
			int start = out.rows.size();
			for (int r = qStart[s]; r < qStart[s + 1]; r++) {
				Hypothesis h = rows.get(r);
				if (candidateOk.test(h)) {
					out.rows.add(h);
					float[] f = x[r];
					if (extra != null) {
						float[] e = extra.apply(s);
						f = Arrays.copyOf(f, f.length + e.length);
						System.arraycopy(e, 0, f, x[r].length, e.length);
					}
					out.x.add(f);
				}
			}
			if (out.rows.size() > start) {
				out.qStart.add(start);
				out.spectrum.add(s);
				out.fold.add(spectrumFold(s));
			}
		}
		return out;
	}

	/**
	 * Stage-wise, exact-first identification ({@code docs/developer/stagewise-exact-first-identification-design.md},
	 * revision 2).
	 * <ol>
	 * <li>Stage 1: each spectrum's exact-mass candidates (closed peptidoforms and open zero family);
	 * winner by the stage-1 model; one TDC.</li>
	 * <li>Blocking: every spectrum whose stage-1 winner scores at or above tau1 (the lowest stage-1
	 * score at PSM q <= fdr), whatever its label, leaves the pool. Decided at PSM level only.</li>
	 * <li>Stage 2: the remaining spectra's retained-peak candidates; winner by the stage-2 model
	 * (with the spectrum's stage-1 outcome as features); Group-walk over KS-merged peak groups.</li>
	 * <li>Peptides: walked in stage 1; stage-2 units whose key passed stage 1 are not walked again.</li>
	 * </ol>
	 * Unexplained candidates are not reported. RESET roles are per decoy peptide in every stage.
	 */
	private void stagewise(List<Hypothesis> rows, int[] qStart, float[][] x, String[] names, int[] monotone,
			File outDir) throws IOException {
		int spectra = qStart.length - 1;

		// ---- Stage 1
		StageQueries s1 = stageQueries(rows, qStart, x, s -> true, h -> h.explanation == Hypothesis.Explanation.EXACT,
				null);
		int[] q1Start = s1.qStartArray();
		double[] u1 = fitQueries(s1.rows, q1Start, s1.xArray(), s1.foldArray(), monotone, names, "stage1_");
		int n1 = s1.spectrum.size();
		Winner[] w1 = new Winner[n1];
		IntStream.range(0, n1).parallel().forEach(k -> w1[k] = winner(s1.rows, q1Start, u1, k));
		double[] q1 = walkStage(w1, null, "stage1_");
		double tau1 = Double.POSITIVE_INFINITY;
		for (int k = 0; k < n1; k++) {
			if (!Double.isNaN(q1[k]) && q1[k] <= cfg.fdr) {
				tau1 = Math.min(tau1, w1[k].score);
			}
		}
		boolean[] blocked = new boolean[spectra];
		double[] stage1Score = new double[spectra];
		boolean[] hadStage1 = new boolean[spectra];
		Arrays.fill(stage1Score, -50);
		int blockedCount = 0;
		for (int k = 0; k < n1; k++) {
			int s = s1.spectrum.get(k);
			hadStage1[s] = true;
			stage1Score[s] = Math.max(-50, w1[k].score);
			if (w1[k].score >= tau1) {
				blocked[s] = true;
				blockedCount++;
			}
		}
		summary.add(String.format(Locale.ROOT, "stage1_threshold_score\t%.4f", tau1));
		summary.add("stage1_spectra_blocked_from_stage2_[all_labels]\t" + blockedCount);

		// ---- Stage 2
		String[] names2 = concat(names, new String[] { "stage1_best_score", "stage1_had_exact" });
		int[] monotone2 = Arrays.copyOf(monotone, monotone.length + 2);
		StageQueries s2 = stageQueries(rows, qStart, x, s -> !blocked[s],
				h -> h.explanation == Hypothesis.Explanation.PEAK,
				s -> new float[] { (float) stage1Score[s], hadStage1[s] ? 1 : 0 });
		int[] q2Start = s2.qStartArray();
		int n2 = s2.spectrum.size();
		Winner[] w2 = new Winner[n2];
		double[] q2 = new double[n2];
		Map<String, Integer> groupOf2 = new HashMap<>();
		int groups2 = 1;
		if (n2 > 0) {
			double[] u2 = fitQueries(s2.rows, q2Start, s2.xArray(), s2.foldArray(), monotone2, names2, "stage2_");
			IntStream.range(0, n2).parallel().forEach(k -> w2[k] = winner(s2.rows, q2Start, u2, k));
			// Peak groups (CONGA merging), from walked winners' scores, labels pooled.
			Map<String, List<Double>> batchScores = new LinkedHashMap<>();
			for (Winner w : w2) {
				if (!isTrainingDecoy(w.hyp)) {
					batchScores.computeIfAbsent(batchKey(w.hyp), b -> new ArrayList<>()).add(w.score);
				}
			}
			List<String> keys = new ArrayList<>(batchScores.keySet());
			List<double[]> arrays = new ArrayList<>();
			for (String b : keys) {
				arrays.add(batchScores.get(b).stream().mapToDouble(Double::doubleValue).toArray());
			}
			int[] merged = GroupWalk.mergeBatches(arrays, -1, 0.01, 2 * cfg.groupWalkK);
			for (int b = 0; b < keys.size(); b++) {
				groupOf2.put(keys.get(b), merged[b]);
				groups2 = Math.max(groups2, merged[b] + 1);
			}
			summary.add("stage2_fdr_groups\t" + groups2 + " groups from " + keys.size() + " peak batches");
			q2 = walkStage(w2, groupOf2, "stage2_");
		}

		// ---- Peptide level, stage by stage
		PeptideLevel p1 = peptideLevel(Arrays.asList(winnerHyps(w1)), winnerScores(w1), q1, new int[n1], 1, "stage1_");
		Map<String, Double> unitQ = new HashMap<>();
		Set<String> passed1 = new HashSet<>();
		for (int i = 0; i < p1.keys.size(); i++) {
			unitQ.put(p1.keys.get(i), p1.q[i]);
			if (!Double.isNaN(p1.q[i]) && p1.q[i] <= cfg.fdr) {
				passed1.add(p1.keys.get(i));
			}
		}
		if (n2 > 0) {
			List<Hypothesis> ph = new ArrayList<>();
			List<double[]> pr = new ArrayList<>();
			for (int k = 0; k < n2; k++) {
				Hypothesis h = w2[k].hyp;
				String key = (h.target ? "T:" : "D:") + Hypothesis.unit(h.peptide);
				Integer g = groupOf2.get(batchKey(h));
				if (!passed1.contains(key)) {
					ph.add(h);
					pr.add(new double[] { w2[k].score, q2[k], g == null ? 0 : g });
				}
			}
			PeptideLevel p2 = peptideLevel(ph, column(pr, 0), column(pr, 1), intColumn(pr, 2), groups2, "stage2_");
			for (int i = 0; i < p2.keys.size(); i++) {
				if (!Double.isNaN(p2.q[i])) {
					unitQ.merge(p2.keys.get(i), p2.q[i], (a, b) -> Double.isNaN(a) ? b : Math.min(a, b));
				}
			}
		}

		// ---- One reported PSM per spectrum: the stage-1 winner of a blocked spectrum, otherwise the
		// stage-2 winner, otherwise the (unaccepted) stage-1 winner.
		Winner[] reported = new Winner[spectra];
		double[] psmQ = new double[spectra];
		Arrays.fill(psmQ, Double.NaN);
		int[] stageOf = new int[spectra];
		for (int k = 0; k < n1; k++) {
			int s = s1.spectrum.get(k);
			reported[s] = w1[k];
			psmQ[s] = q1[k];
			stageOf[s] = 1;
		}
		for (int k = 0; k < n2; k++) {
			int s = s2.spectrum.get(k);
			reported[s] = w2[k];
			psmQ[s] = q2[k];
			stageOf[s] = 2;
		}
		int[][] kept = new int[3][2];
		int[][] accepted = new int[3][2];
		Set<String>[] peptides = new Set[] { new HashSet<String>(), new HashSet<String>(), new HashSet<String>() };
		for (int s = 0; s < spectra; s++) {
			Winner w = reported[s];
			if (w == null || isTrainingDecoy(w.hyp) || Double.isNaN(psmQ[s]) || psmQ[s] > cfg.fdr) {
				continue;
			}
			int side = w.hyp.target ? 0 : 1;
			accepted[stageOf[s]][side]++;
			String key = (w.hyp.target ? "T:" : "D:") + Hypothesis.unit(w.hyp.peptide);
			Double pq = unitQ.get(key);
			if (pq != null && !Double.isNaN(pq) && pq <= cfg.fdr) {
				kept[stageOf[s]][side]++;
				peptides[stageOf[s]].add(key);
			}
		}
		Set<String> allPeptides = new HashSet<>(peptides[1]);
		allPeptides.addAll(peptides[2]);
		int pepT = 0;
		int pepD = 0;
		for (String k : allPeptides) {
			if (k.startsWith("T:")) {
				pepT++;
			} else {
				pepD++;
			}
		}
		for (int st = 1; st <= 2; st++) {
			summary.add("stage" + st + "_psms_at_fdr_[targets,estimating_decoys]\t" + Arrays.toString(accepted[st]));
			summary.add("stage" + st + "_psms_kept_[targets,estimating_decoys]\t" + Arrays.toString(kept[st]));
		}
		int keptT = kept[1][0] + kept[2][0];
		int keptD = kept[1][1] + kept[2][1];
		summary.add("fdr\t" + cfg.fdr + " per stage (stage 1 TDC, stage 2 Group-walk K=" + cfg.groupWalkK
				+ "), estimate (" + decoyWeight() + "*D+1)/T, RESET " + cfg.reset + "; union not guaranteed (design 3, T3)");
		summary.add("identifications_psms_kept\t[" + keptT + ", " + keptD + "]");
		summary.add("identifications_peptides\t[" + pepT + ", " + pepD + "]");
		summary.add(String.format(Locale.ROOT, "estimated_psm_fdr_union\t%.5f", decoyWeight() * keptD / Math.max(1, keptT)));
		summary.add(String.format(Locale.ROOT, "estimated_peptide_fdr_union\t%.5f", decoyWeight() * pepD / Math.max(1, pepT)));
		for (String line : summary) {
			log.accept("Unified: " + line.replace('\t', ' '));
		}
		outDir.mkdirs();
		batches.writeReport(new File(outDir, "unified_peaks.tsv"));
		try (PrintWriter w = new PrintWriter(new File(outDir, "unified_summary.tsv"), "UTF-8")) {
			w.println("name\tvalue");
			for (String line : summary) {
				w.println(line);
			}
		}
		writeWinners(reported, psmQ, unitQ, groupOf2, outDir);
	}

	/**
	 * PSM q-values of one stage's winners: TDC (one group) or Group-walk over {@code groupOf}, with
	 * estimating decoys only. NaN for training decoys.
	 */
	private double[] walkStage(Winner[] w, Map<String, Integer> groupOf, String label) {
		List<Integer> items = new ArrayList<>();
		for (int k = 0; k < w.length; k++) {
			if (!isTrainingDecoy(w[k].hyp)) {
				items.add(k);
			}
		}
		int m = items.size();
		double[] sc = new double[m];
		boolean[] tg = new boolean[m];
		int[] gr = new int[m];
		int groups = 1;
		int[] counts = new int[3];
		for (int i = 0; i < m; i++) {
			Winner x = w[items.get(i)];
			sc[i] = x.score;
			tg[i] = x.hyp.target;
			Integer g = groupOf == null ? null : groupOf.get(batchKey(x.hyp));
			gr[i] = g == null ? 0 : g;
			groups = Math.max(groups, gr[i] + 1);
			counts[x.hyp.target ? 0 : 1]++;
		}
		counts[2] = w.length - m;
		double[] q = GroupWalk.qValues(sc, tg, gr, groups, cfg.groupWalkK, decoyWeight(), 1.0);
		double[] out = new double[w.length];
		Arrays.fill(out, Double.NaN);
		for (int i = 0; i < m; i++) {
			out[items.get(i)] = q[i];
		}
		summary.add(label + "winners_[targets,estimating_decoys,training_decoys]\t" + Arrays.toString(counts));
		return out;
	}

	private static Hypothesis[] winnerHyps(Winner[] w) {
		Hypothesis[] out = new Hypothesis[w.length];
		for (int i = 0; i < w.length; i++) {
			out[i] = w[i].hyp;
		}
		return out;
	}

	private static double[] winnerScores(Winner[] w) {
		double[] out = new double[w.length];
		for (int i = 0; i < w.length; i++) {
			out[i] = w[i].score;
		}
		return out;
	}

	// ------------------------------------------------------------ Group-walk and outputs

	private void walkAndWrite(Winner[] winners, File outDir) throws IOException {
		int spectra = winners.length;
		// Batch of each winner; training decoys are out (they trained the model).
		Map<String, List<Double>> batchScores = new LinkedHashMap<>();
		for (Winner w : winners) {
			if (w == null || isTrainingDecoy(w.hyp)) {
				continue;
			}
			if (w.hyp.explanation == Hypothesis.Explanation.UNEXPLAINED && !cfg.reportUnexplained) {
				continue;
			}
			batchScores.computeIfAbsent(batchKey(w.hyp), k -> new ArrayList<>()).add(w.score);
		}
		// Groups (CONGA): the exact batch fixed; the other batches merged when a KS test cannot tell
		// their winning-score distributions apart, and every group kept at >= 2K winners. Scores of
		// targets and decoys pooled, so no label is read.
		List<String> batchKeys = new ArrayList<>(batchScores.keySet());
		List<double[]> scoreArrays = new ArrayList<>();
		for (String k : batchKeys) {
			scoreArrays.add(batchScores.get(k).stream().mapToDouble(Double::doubleValue).toArray());
		}
		int[] merged = GroupWalk.mergeBatches(scoreArrays, batchKeys.indexOf("EXACT"), 0.01, 2 * cfg.groupWalkK);
		Map<String, Integer> groupOf = new LinkedHashMap<>();
		int groupCount = 0;
		for (int g : merged) {
			groupCount = Math.max(groupCount, g + 1);
		}
		List<String> groupNames = new ArrayList<>();
		for (int g = 0; g < groupCount; g++) {
			List<String> members = new ArrayList<>();
			for (int b = 0; b < batchKeys.size(); b++) {
				if (merged[b] == g) {
					members.add(batchKeys.get(b));
				}
			}
			groupNames.add(members.size() == 1 ? members.get(0) : "GROUP_" + g + "(" + members.size() + "_batches)");
		}
		for (int b = 0; b < batchKeys.size(); b++) {
			groupOf.put(batchKeys.get(b), merged[b]);
		}
		summary.add("fdr_groups\t" + groupCount + " groups from " + batchKeys.size() + " batches: "
				+ String.join(", ", groupNames));
		List<Integer> items = new ArrayList<>();
		for (int s = 0; s < spectra; s++) {
			Winner w = winners[s];
			if (w != null && !isTrainingDecoy(w.hyp) && groupOf.containsKey(batchKey(w.hyp))) {
				items.add(s);
			}
		}
		int m = items.size();
		double[] sc = new double[m];
		boolean[] tg = new boolean[m];
		int[] gr = new int[m];
		for (int i = 0; i < m; i++) {
			Winner w = winners[items.get(i)];
			sc[i] = w.score;
			tg[i] = w.hyp.target;
			gr[i] = groupOf.get(batchKey(w.hyp));
		}
		double[] q = GroupWalk.qValues(sc, tg, gr, groupNames.size(), cfg.groupWalkK, decoyWeight(), 1.0);
		double[] psmQ = new double[spectra];
		Arrays.fill(psmQ, Double.NaN);
		for (int i = 0; i < m; i++) {
			psmQ[items.get(i)] = q[i];
		}

		// Peptide level: one unit per peptide from all its walked PSM winners, scored by the
		// peptide-level model (support evidence), then walked per batch.
		List<Hypothesis> unitPsms = new ArrayList<>();
		List<double[]> unitRows = new ArrayList<>(); // {score, q, group}
		for (int i = 0; i < m; i++) {
			unitPsms.add(winners[items.get(i)].hyp);
			unitRows.add(new double[] { sc[i], q[i], gr[i] });
		}
		for (Winner w : winners) {
			if (w != null && isTrainingDecoy(w.hyp) && groupOf.containsKey(batchKey(w.hyp))) {
				unitPsms.add(w.hyp);
				unitRows.add(new double[] { w.score, Double.NaN, groupOf.get(batchKey(w.hyp)) });
			}
		}
		PeptideLevel pl = peptideLevel(unitPsms, column(unitRows, 0), column(unitRows, 1), intColumn(unitRows, 2),
				groupNames.size(), "");
		List<String> unitKeys = pl.keys;
		double[] uq = pl.q;
		int[] ug = pl.group;
		boolean[] ut = pl.target;
		Map<String, Double> unitQ = new HashMap<>();
		for (int i = 0; i < unitKeys.size(); i++) {
			unitQ.put(unitKeys.get(i), uq[i]);
		}

		// Summary per group.
		int g = groupNames.size();
		int[][] winnersBy = new int[g][3];
		int[][] psmsBy = new int[g][2];
		int[][] keptBy = new int[g][2];
		for (int i = 0; i < m; i++) {
			Hypothesis h = winners[items.get(i)].hyp;
			int side = h.target ? 0 : 1;
			winnersBy[gr[i]][side]++;
			if (q[i] <= cfg.fdr) {
				psmsBy[gr[i]][side]++;
				Double pq = unitQ.get((h.target ? "T:" : "D:") + Hypothesis.unit(h.peptide));
				if (pq != null && pq <= cfg.fdr) {
					keptBy[gr[i]][side]++;
				}
			}
		}
		for (Winner w : winners) {
			if (w != null && isTrainingDecoy(w.hyp) && groupOf.containsKey(batchKey(w.hyp))) {
				winnersBy[groupOf.get(batchKey(w.hyp))][2]++;
			}
		}
		int[] pepGroup = new int[g * 2];
		for (int i = 0; i < unitKeys.size(); i++) {
			if (uq[i] <= cfg.fdr) {
				pepGroup[ug[i] * 2 + (ut[i] ? 0 : 1)]++;
			}
		}
		summary.add("fdr\t" + cfg.fdr + " (Group-walk, K=" + cfg.groupWalkK + ", estimate (" + decoyWeight()
				+ "*D+1)/T, RESET " + cfg.reset + ")");
		summary.add("mode\tclosed_only=" + cfg.closedOnly + " engine_features=" + cfg.engineFeatures
				+ " max_training_queries=" + cfg.maxTrainingQueries);
		summary.add("unexplained_reported\t" + cfg.reportUnexplained);
		int[] totPsm = new int[2];
		int[] totKept = new int[2];
		int[] totPep = new int[2];
		for (int k = 0; k < g; k++) {
			String name = groupNames.get(k).toLowerCase(Locale.ROOT);
			summary.add(name + "_winners_[targets,estimating_decoys,training_decoys]\t" + Arrays.toString(winnersBy[k]));
			summary.add(name + "_psms_at_fdr\t" + Arrays.toString(psmsBy[k]));
			summary.add(name + "_psms_kept\t" + Arrays.toString(keptBy[k]));
			summary.add(name + "_peptides_at_fdr\t[" + pepGroup[k * 2] + ", " + pepGroup[k * 2 + 1] + "]");
			for (int side = 0; side < 2; side++) {
				totPsm[side] += psmsBy[k][side];
				totKept[side] += keptBy[k][side];
				totPep[side] += pepGroup[k * 2 + side];
			}
		}
		summary.add("identifications_psms_kept\t" + Arrays.toString(totKept));
		summary.add("identifications_peptides\t" + Arrays.toString(totPep));
		summary.add(String.format(Locale.ROOT, "estimated_psm_fdr\t%.5f",
				(decoyWeight() * totKept[1]) / Math.max(1, totKept[0])));
		summary.add(String.format(Locale.ROOT, "estimated_peptide_fdr\t%.5f",
				(decoyWeight() * totPep[1]) / Math.max(1, totPep[0])));
		for (String line : summary) {
			log.accept("Unified: " + line.replace('\t', ' '));
		}

		outDir.mkdirs();
		batches.writeReport(new File(outDir, "unified_peaks.tsv"));
		try (PrintWriter w = new PrintWriter(new File(outDir, "unified_summary.tsv"), "UTF-8")) {
			w.println("name\tvalue");
			for (String line : summary) {
				w.println(line);
			}
		}
		writeWinners(winners, psmQ, unitQ, groupOf, outDir);
	}

	/**
	 * Peptide-level evidence, computed from a peptide's walked PSM winners identically for target and
	 * decoy peptides. Recurrence is evidence that a peptide is present, which is the peptide-level
	 * question; an absent target peptide that a random match hits does not recur, like a decoy. The
	 * PSM score stays spectrum-intrinsic, so PSM-level FDR is not lifted by recurrence.
	 */
	static final String[] PEPTIDE_FEATURES = { "best_psm_score", "second_psm_score", "top3_mean_psm_score",
			"log_n_psms", "log_n_psms_at_fdr", "run_fraction", "log_n_runs", "n_charges", "peptide_length",
			"missed_cleavages", "found_by_both_fraction", "best_is_exact", "best_is_peak", "best_is_unexplained",
			"n_explanations" };

	private static final class PeptideLevel {
		List<String> keys;
		double[] q;
		int[] group;
		boolean[] target;
	}

	/**
	 * Peptide units ({@code T:}/{@code D:} + base sequence, I = L) from a stage's PSM winners, scored
	 * by a semi-supervised, cross-validated model over {@link #PEPTIDE_FEATURES}, then Group-walk over
	 * the batches of each unit's best PSM.
	 * <p>
	 * RESET at the peptide level (review R3): decoy roles are assigned per decoy peptide, so a decoy
	 * unit is wholly training or wholly estimating. Training-decoy units are the model's negatives and
	 * never enter the walk; estimating-decoy units never train and are the walk's only decoys.
	 * Positives are pseudo-target units (targets and estimating decoys) at the seeding threshold,
	 * which uses about 3 false pseudo-targets per training decoy. Without RESET every decoy unit
	 * trains and estimates (Percolator style).
	 *
	 * @param psms the stage's winners, training decoys included
	 * @param q    PSM q-values; NaN for training decoys (they are not walked)
	 */
	private PeptideLevel peptideLevel(List<Hypothesis> psms, double[] sc, double[] q, int[] gr, int groups,
			String label) {
		Map<String, List<Integer>> byUnit = new HashMap<>();
		Set<String> allRuns = new HashSet<>();
		for (int i = 0; i < psms.size(); i++) {
			Hypothesis h = psms.get(i);
			byUnit.computeIfAbsent((h.target ? "T:" : "D:") + Hypothesis.unit(h.peptide), k -> new ArrayList<>()).add(i);
			allRuns.add(specFile.get(h.spectrum));
		}
		// "PSMs at fdr" as a score threshold, so a training-decoy unit (never walked) is counted like
		// any other: the lowest score among walked PSMs at q <= fdr.
		double tau = Double.POSITIVE_INFINITY;
		for (int i = 0; i < psms.size(); i++) {
			if (!Double.isNaN(q[i]) && q[i] <= cfg.fdr) {
				tau = Math.min(tau, sc[i]);
			}
		}
		List<String> keys = new ArrayList<>(byUnit.keySet());
		java.util.Collections.sort(keys);
		int n = keys.size();
		double minScore = Double.POSITIVE_INFINITY;
		for (double s : sc) {
			minScore = Math.min(minScore, s);
		}
		float[][] x = new float[n][];
		double[] best = new double[n];
		int[] group = new int[n];
		boolean[] target = new boolean[n];
		boolean[] trainingUnit = new boolean[n];
		for (int u = 0; u < n; u++) {
			List<Integer> idx = byUnit.get(keys.get(u));
			idx.sort((a, b) -> Double.compare(sc[b], sc[a]));
			Hypothesis top = psms.get(idx.get(0));
			trainingUnit[u] = !top.target && top.resetTraining;
			Set<String> runs = new HashSet<>();
			Set<Integer> charges = new HashSet<>();
			Set<Hypothesis.Explanation> explanations = new HashSet<>();
			int atFdr = 0;
			int both = 0;
			double top3 = 0;
			for (int k = 0; k < idx.size(); k++) {
				int i = idx.get(k);
				Hypothesis h = psms.get(i);
				runs.add(specFile.get(h.spectrum));
				charges.add(h.charge);
				explanations.add(h.explanation);
				if (sc[i] >= tau) {
					atFdr++;
				}
				if (h.foundByClosed && h.foundByOpen) {
					both++;
				}
				if (k < 3) {
					top3 += sc[i];
				}
			}
			best[u] = sc[idx.get(0)];
			group[u] = gr[idx.get(0)];
			target[u] = keys.get(u).startsWith("T:");
			x[u] = new float[] { (float) best[u], (float) (idx.size() > 1 ? sc[idx.get(1)] : minScore),
					(float) (top3 / Math.min(3, idx.size())), (float) Math.log1p(idx.size()), (float) Math.log1p(atFdr),
					runs.size() / (float) Math.max(1, allRuns.size()), (float) Math.log1p(runs.size()), charges.size(),
					top.peptide.length(), Math.max(0, top.missedCleavages), both / (float) idx.size(),
					top.explanation == Hypothesis.Explanation.EXACT ? 1 : 0,
					top.explanation == Hypothesis.Explanation.PEAK ? 1 : 0,
					top.explanation == Hypothesis.Explanation.UNEXPLAINED ? 1 : 0, explanations.size() };
		}
		Set<String> up = new HashSet<>(Arrays.asList("best_psm_score", "second_psm_score", "top3_mean_psm_score",
				"log_n_psms_at_fdr", "run_fraction", "found_by_both_fraction"));
		int[] monotone = new int[PEPTIDE_FEATURES.length];
		for (int k = 0; k < monotone.length; k++) {
			monotone[k] = up.contains(PEPTIDE_FEATURES[k]) ? 1 : 0;
		}
		int[] fold = new int[n];
		for (int u = 0; u < n; u++) {
			fold[u] = (int) Math.floorMod(stableHash("pep:" + keys.get(u).substring(2)), (long) cfg.folds);
		}
		double[] score = new double[n];
		for (int f = 0; f < cfg.folds; f++) {
			double[] working = best.clone();
			ListwiseGbdt.Model model = null;
			for (int it = 0; it < cfg.iterations; it++) {
				// Seeds over the training folds. Under RESET: negatives are training-decoy units,
				// positives pseudo-target units (targets and estimating decoys, indistinguishable to the
				// trainer) above the threshold where ~3 false pseudo-targets per training decoy reach
				// trainQ. Without RESET: target units vs every decoy unit, decoys / targets.
				List<Integer> train = new ArrayList<>();
				for (int u = 0; u < n; u++) {
					if (fold[u] != f) {
						train.add(u);
					}
				}
				final double[] wk = working;
				train.sort((a, b) -> Double.compare(wk[b], wk[a]));
				double[] fdr = new double[train.size()];
				int pseudo = 0;
				int negative = 0;
				for (int k = 0; k < train.size(); k++) {
					int u = train.get(k);
					boolean isNegative = cfg.reset ? trainingUnit[u] : !target[u];
					if (isNegative) {
						negative++;
					} else {
						pseudo++;
					}
					fdr[k] = pseudo == 0 ? 1.0 : Math.min(1.0, (cfg.reset ? 3.0 : 1.0) * negative / pseudo);
				}
				for (int k = fdr.length - 2; k >= 0; k--) {
					fdr[k] = Math.min(fdr[k], fdr[k + 1]);
				}
				List<Integer> rows = new ArrayList<>();
				List<Integer> positive = new ArrayList<>();
				for (int k = 0; k < train.size(); k++) {
					int u = train.get(k);
					boolean isNegative = cfg.reset ? trainingUnit[u] : !target[u];
					if (isNegative) {
						rows.add(u);
						positive.add(-1);
					} else if (fdr[k] <= cfg.trainQ) {
						rows.add(u);
						positive.add(1);
					}
				}
				long pos = positive.stream().filter(v -> v == 1).count();
				if (pos < 50 || rows.size() - pos < 50) {
					break;
				}
				float[][] tx = new float[rows.size()][];
				int[] qs = new int[rows.size() + 1];
				int[] tp = new int[rows.size()];
				for (int r = 0; r < rows.size(); r++) {
					tx[r] = x[rows.get(r)];
					qs[r] = r;
					tp[r] = positive.get(r) == 1 ? r : -1;
				}
				qs[rows.size()] = rows.size();
				model = ListwiseGbdt.train(tx, qs, tp, monotone, cfg.gbdt);
				if (it == cfg.iterations - 1) {
					summary.add(String.format(Locale.ROOT, "%speptide_fold_%d_[positive,negative]_units\t[%d, %d]", label,
							f + 1, pos, rows.size() - pos));
				}
				for (int u : train) {
					working[u] = model.utility(x[u]);
				}
			}
			for (int u = 0; u < n; u++) {
				if (fold[u] == f) {
					// Trees give one value per leaf, so many units tie; the best PSM score orders them.
					score[u] = model == null ? best[u] : model.utility(x[u]) + 1e-6 * best[u];
				}
			}
		}
		// The walks: target units and estimating-decoy units only (training-decoy units trained).
		List<Integer> walked = new ArrayList<>();
		for (int u = 0; u < n; u++) {
			if (!(cfg.reset && trainingUnit[u])) {
				walked.add(u);
			}
		}
		int m = walked.size();
		double[] ws = new double[m];
		double[] wb = new double[m];
		boolean[] wt = new boolean[m];
		int[] wg = new int[m];
		List<double[]> a1 = new ArrayList<>();
		List<double[]> a2 = new ArrayList<>();
		for (int k = 0; k < m; k++) {
			int u = walked.get(k);
			ws[k] = score[u];
			wb[k] = best[u];
			wt[k] = target[u];
			wg[k] = group[u];
			a1.add(new double[] { score[u], target[u] ? 1 : 0 });
			a2.add(new double[] { best[u], target[u] ? 1 : 0 });
		}
		double[] wModelQ = GroupWalk.qValues(ws, wt, wg, groups, cfg.groupWalkK, decoyWeight(), 1.0);
		double[] wBestQ = GroupWalk.qValues(wb, wt, wg, groups, cfg.groupWalkK, decoyWeight(), 1.0);
		double[] modelQ = new double[n];
		double[] bestQ = new double[n];
		Arrays.fill(modelQ, Double.NaN);
		Arrays.fill(bestQ, Double.NaN);
		int bestPass = 0;
		int modelPass = 0;
		for (int k = 0; k < m; k++) {
			int u = walked.get(k);
			modelQ[u] = wModelQ[k];
			bestQ[u] = wBestQ[k];
			if (target[u] && wBestQ[k] <= cfg.fdr) {
				bestPass++;
			}
			if (target[u] && wModelQ[k] <= cfg.fdr) {
				modelPass++;
			}
		}
		summary.add(String.format(Locale.ROOT, "%speptide_units_[targets,estimating_decoys,training_decoys]\t[%d, %d, %d]",
				label, countTrue(target), m - countTrue(target), n - m));
		summary.add(String.format(Locale.ROOT, "%speptide_auc_[model,best_psm_only]\t[%.4f, %.4f]", label, auc(a1),
				auc(a2)));
		summary.add(label + "peptides_at_fdr_[model,best_psm_only]\t[" + modelPass + ", " + bestPass + "]");
		if (unitTable != null) {
			File table = new File(unitTable.getParentFile(), label + unitTable.getName());
			try (PrintWriter w = new PrintWriter(new BufferedWriter(new OutputStreamWriter(
					new GZIPOutputStream(new FileOutputStream(table)), StandardCharsets.UTF_8)))) {
				w.println("unit\tgroup\tfold\tmodel_score\tmodel_q\tbest_psm_q\t" + String.join("\t", PEPTIDE_FEATURES));
				for (int u = 0; u < n; u++) {
					StringBuilder sb = new StringBuilder(keys.get(u)).append('\t').append(group[u]).append('\t')
							.append(fold[u]).append('\t').append(String.format(Locale.ROOT, "%.5f", score[u]))
							.append('\t').append(String.format(Locale.ROOT, "%.6f", modelQ[u])).append('\t')
							.append(String.format(Locale.ROOT, "%.6f", bestQ[u]));
					for (float v : x[u]) {
						sb.append('\t').append(v);
					}
					w.println(sb);
				}
			} catch (IOException e) {
				log.accept("Unified: WARNING could not write the peptide unit table: " + e.getMessage());
			}
		}
		PeptideLevel out = new PeptideLevel();
		out.keys = keys;
		out.q = modelQ;
		out.group = group;
		out.target = target;
		return out;
	}

	private static double[] column(List<double[]> rows, int k) {
		double[] out = new double[rows.size()];
		for (int i = 0; i < out.length; i++) {
			out[i] = rows.get(i)[k];
		}
		return out;
	}

	private static int[] intColumn(List<double[]> rows, int k) {
		int[] out = new int[rows.size()];
		for (int i = 0; i < out.length; i++) {
			out[i] = (int) rows.get(i)[k];
		}
		return out;
	}

	private static int countTrue(boolean[] v) {
		int c = 0;
		for (boolean b : v) {
			if (b) {
				c++;
			}
		}
		return c;
	}

	/** Excluded from the walks: a decoy that trained the model under RESET. Without RESET none is. */
	private boolean isTrainingDecoy(Hypothesis h) {
		return cfg.reset && !h.target && h.resetTraining;
	}

	private double decoyWeight() {
		return cfg.reset ? 2.0 : 1.0;
	}

	private void writeWinners(Winner[] winners, double[] psmQ, Map<String, Double> unitQ, Map<String, Integer> groupOf,
			File outDir) throws IOException {
		Map<String, String> proteins = new LinkedHashMap<>();
		try (PrintWriter psms = new PrintWriter(new BufferedWriter(new OutputStreamWriter(
				new GZIPOutputStream(new FileOutputStream(new File(outDir, "unified_psms.tsv.gz"))),
				StandardCharsets.UTF_8), 1 << 20));
				PrintWriter arms = new PrintWriter(new BufferedWriter(new OutputStreamWriter(
						new GZIPOutputStream(new FileOutputStream(new File(outDir, "unified_arms.tsv.gz"))),
						StandardCharsets.UTF_8), 1 << 20))) {
			psms.println("file\tscan\tpeptide\tis_target\treset_role\texplanation\tbatch\tdelta\tisotope_offset\t"
					+ "found_by_closed\tfound_by_open\tutility\tscore\tpsm_q\tpeptide_q\tkept\tproteins");
			arms.println("arm\tfile\tscan\tpeptide\tis_target\tclass\tq_value\taccepted");
			for (int s = 0; s < winners.length; s++) {
				Winner w = winners[s];
				if (w == null) {
					continue;
				}
				Hypothesis h = w.hyp;
				String role = h.target ? "target" : isTrainingDecoy(h) ? "training_decoy" : "estimating_decoy";
				Double pq = Double.isNaN(psmQ[s]) ? null
						: unitQ.get((h.target ? "T:" : "D:") + Hypothesis.unit(h.peptide));
				boolean kept = h.target && !Double.isNaN(psmQ[s]) && psmQ[s] <= cfg.fdr && pq != null && pq <= cfg.fdr;
				String prot = String.join(";", Arrays.asList(h.proteins).stream()
						.filter(p -> p != null && !p.isBlank()).toArray(String[]::new));
				psms.println(specFile.get(h.spectrum) + "\t" + specScan.get(h.spectrum) + "\t" + h.peptide + "\t"
						+ h.target + "\t" + role + "\t" + h.explanation + "\t" + batchKey(h) + "\t"
						+ String.format(Locale.ROOT, "%.5f", h.foundByClosed ? h.closedMassDiff : h.fragmentDelta) + "\t"
						+ h.isotopeOffset + "\t" + h.foundByClosed + "\t" + h.foundByOpen + "\t"
						+ String.format(Locale.ROOT, "%.4f\t%.4f", w.utility, w.score) + "\t"
						+ (Double.isNaN(psmQ[s]) ? "" : String.format(Locale.ROOT, "%.6f", psmQ[s])) + "\t"
						+ (pq == null ? "" : String.format(Locale.ROOT, "%.6f", pq)) + "\t" + kept + "\t" + prot);
				if (!Double.isNaN(psmQ[s]) && psmQ[s] <= 0.05 && !"training_decoy".equals(role)) {
					String cls = h.explanation == Hypothesis.Explanation.EXACT ? "exact"
							: h.explanation == Hypothesis.Explanation.PEAK ? "open_peak" : "open_unexplained";
					// Decoys are marked as targets would be, so decoy ratios can be read (x2: RESET).
					boolean accepted = psmQ[s] <= cfg.fdr && pq != null && pq <= cfg.fdr;
					arms.println((cfg.stagewise ? "stagewise" : "unified") + "\t" + specFile.get(h.spectrum) + "\t" + specScan.get(h.spectrum) + "\t"
							+ h.peptide + "\t" + h.target + "\t" + cls + "\t" + psmQ[s] + "\t" + accepted);
					proteins.putIfAbsent(h.peptide, prot);
				}
			}
		}
		try (PrintWriter w = new PrintWriter(new BufferedWriter(new OutputStreamWriter(
				new GZIPOutputStream(new FileOutputStream(new File(outDir, "unified_arm_peptides.tsv.gz"))),
				StandardCharsets.UTF_8)))) {
			w.println("peptide\tproteins");
			for (Map.Entry<String, String> e : proteins.entrySet()) {
				w.println(e.getKey() + "\t" + e.getValue());
			}
		}
	}

	// ------------------------------------------------------------ helpers

	private static int[] monotone(String[] names) {
		Set<String> up = new HashSet<>(Arrays.asList("u_matched", "u_int_tic", "u_int_top50", "e_matched", "e_int_tic",
				"e_hyperscore", "u_hyperscore", "p_entropy", "p_pearson", "p_observed"));
		Set<String> down = new HashSet<>(Arrays.asList("expl_error_sigma", "p_rt_error"));
		int[] m = new int[names.length];
		for (int i = 0; i < names.length; i++) {
			m[i] = up.contains(names[i]) ? 1 : down.contains(names[i]) ? -1 : 0;
		}
		return m;
	}

	private static String[] concat(String[] a, String[] b) {
		String[] out = Arrays.copyOf(a, a.length + b.length);
		System.arraycopy(b, 0, out, a.length, b.length);
		return out;
	}

	static double median(double[] v) {
		if (v.length == 0) {
			return 0;
		}
		double[] c = v.clone();
		Arrays.sort(c);
		return c.length % 2 == 1 ? c[c.length / 2] : (c[c.length / 2 - 1] + c[c.length / 2]) / 2;
	}

	/** FNV-1a, deterministic across JVMs. */
	static long stableHash(String value) {
		long hash = 0xcbf29ce484222325L;
		for (int i = 0; i < value.length(); i++) {
			hash ^= value.charAt(i);
			hash *= 0x100000001b3L;
		}
		return hash;
	}
}
