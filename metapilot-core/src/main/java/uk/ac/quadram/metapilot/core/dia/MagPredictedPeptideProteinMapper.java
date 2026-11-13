package uk.ac.quadram.metapilot.core.dia;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.concurrent.ConcurrentMap;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import uk.ac.quadram.metapilot.dbSearch.diann.DiaNNPrecursor;
import uk.ac.quadram.metapilot.task.io.pep.MetaPSM;

/**
 * Maps identified peptide sequences to MAG genomes/proteins using predicted TSVs.
 *
 * <p>Predicted TSV format: column 0 = protein ID, column 1 = peptide sequence.
 */
public final class MagPredictedPeptideProteinMapper {

	private static final Logger LOGGER = LogManager.getLogger(MagPredictedPeptideProteinMapper.class);
	private static final String PREDICTED_TSV_SUFFIX = ".predicted.tsv";

	private MagPredictedPeptideProteinMapper() {
	}

	public static final class Result {
		public final ConcurrentHashMap<String, HashMap<String, HashSet<String>>> genomeProteinPeptides;
		public final HashSet<String> matchedPeptides;

		private Result(ConcurrentHashMap<String, HashMap<String, HashSet<String>>> genomeProteinPeptides,
				HashSet<String> matchedPeptides) {
			this.genomeProteinPeptides = genomeProteinPeptides;
			this.matchedPeptides = matchedPeptides;
		}
	}

	public static Result mapPeptides(Set<String> peptides, File predictedFolder, Set<String> allowedGenomes,
			int threadCount) {
		return mapPeptides(peptides, predictedFolder, allowedGenomes, null, threadCount);
	}

	public static Result mapPeptides(Set<String> peptides, File predictedFolder, Set<String> allowedGenomes,
			Set<String> allowedProteins, int threadCount) {
		ConcurrentHashMap<String, HashMap<String, HashSet<String>>> genomeProteinPeptides = new ConcurrentHashMap<>();
		if (peptides == null || peptides.isEmpty()) {
			return new Result(genomeProteinPeptides, new HashSet<>());
		}
		if (predictedFolder == null || !predictedFolder.isDirectory()) {
			LOGGER.error("Predicted peptide folder not found: {}",
					predictedFolder == null ? "null" : predictedFolder.getAbsolutePath());
			return new Result(genomeProteinPeptides, new HashSet<>());
		}

		File[] files = predictedFolder.listFiles(file -> file.isFile() && file.getName().endsWith(PREDICTED_TSV_SUFFIX)
				&& isAllowedGenome(file, allowedGenomes));
		if (files == null || files.length == 0) {
			LOGGER.error("No predicted peptide files found in: {}", predictedFolder.getAbsolutePath());
			return new Result(genomeProteinPeptides, new HashSet<>());
		}

		Set<String> peptideSnapshot = Collections.unmodifiableSet(new HashSet<>(peptides));
		Set<String> matchedPeptides = ConcurrentHashMap.newKeySet();
		ExecutorService executor = Executors.newFixedThreadPool(Math.max(1, Math.min(threadCount, files.length)));
		for (File file : files) {
			executor.execute(() -> mapFile(file, peptideSnapshot, allowedProteins, matchedPeptides, genomeProteinPeptides));
		}
		executor.shutdown();

		try {
			if (!executor.awaitTermination(1, TimeUnit.HOURS)) {
				LOGGER.error("Peptide-to-genome mapping timed out after 1 hour while reading {}", predictedFolder);
				return null;
			}
		} catch (InterruptedException e) {
			Thread.currentThread().interrupt();
			LOGGER.error("Peptide-to-genome mapping was interrupted while reading {}", predictedFolder, e);
			return null;
		}

		return new Result(genomeProteinPeptides, new HashSet<>(matchedPeptides));
	}

	public static ConcurrentHashMap<String, Set<String>> mapPeptideProteins(Set<String> peptides, File predictedFolder,
			Set<String> allowedGenomes, int threadCount) {
		ConcurrentHashMap<String, Set<String>> peptideProteins = new ConcurrentHashMap<>();
		if (peptides == null || peptides.isEmpty()) {
			return peptideProteins;
		}
		for (String peptide : peptides) {
			if (peptide != null && !peptide.isEmpty()) {
				peptideProteins.put(peptide, ConcurrentHashMap.newKeySet());
			}
		}
		if (peptideProteins.isEmpty() || predictedFolder == null || !predictedFolder.isDirectory()) {
			return peptideProteins;
		}

		File[] files = predictedFolder.listFiles();
		if (files != null && files.length > 0) {
			ExecutorService executor = Executors.newFixedThreadPool(Math.min(Math.max(1, threadCount), files.length));
			for (File file : files) {
				if (!file.isFile() || !file.getName().endsWith(PREDICTED_TSV_SUFFIX)) {
					continue;
				}
				if (!isAllowedGenome(file, allowedGenomes)) {
					continue;
				}
				executor.execute(() -> mapProteinAssignments(file, peptideProteins));
			}
			executor.shutdown();
			try {
				executor.awaitTermination(1, TimeUnit.HOURS);
			} catch (InterruptedException e) {
				Thread.currentThread().interrupt();
				LOGGER.warn("Interrupted while reading predicted TSVs");
			}
		}
		return peptideProteins;
	}

    /** Final-export mapping retains the original decoy background and unmapped targets
     * (with empty associations) for diagnostic exclusion. Upstream selection mapping is unchanged. */
    public static DiaNNPrecursor[] mapForFinalisation(DiaNNPrecursor[] precursors, File predictedFolder,
            Set<String> allowedGenomes) {
        if (predictedFolder == null || !predictedFolder.isDirectory())
            throw new IllegalStateException("Seed finalisation requires predicted mapping resources");
        if (precursors == null) return new DiaNNPrecursor[0];
        DiaNNPrecursor[] targets = Arrays.stream(precursors).filter(p -> p != null && !p.getIsDecoy())
                .toArray(DiaNNPrecursor[]::new);
        List<DiaNNPrecursor> mapped = mapPrecursorsToProteins(targets, predictedFolder, allowedGenomes);
        if (mapped == null) throw new IllegalStateException("Seed finalisation mapping failed");
        Set<DiaNNPrecursor> matched = Collections.newSetFromMap(new java.util.IdentityHashMap<>());
        matched.addAll(mapped);
        for (DiaNNPrecursor p : targets) if (!matched.contains(p)) p.setProGroups(new String[0]);
        return precursors;
    }

	public static List<DiaNNPrecursor> mapPrecursorsToProteins(DiaNNPrecursor[] precursors, File predictedFolder) {
		return mapPrecursorsToProteins(precursors, predictedFolder, null);
	}

	public static List<DiaNNPrecursor> mapPrecursorsToProteins(DiaNNPrecursor[] precursors, File predictedFolder,
			Set<String> allowedGenomes) {
		if (precursors == null || precursors.length == 0 || predictedFolder == null || !predictedFolder.isDirectory()) {
			return precursors != null && precursors.length > 0
					? new ArrayList<>(Arrays.asList(precursors))
					: new ArrayList<>();
		}
		HashSet<String> peptides = new HashSet<>();
		for (DiaNNPrecursor precursor : precursors) {
			String sequence = precursor != null ? precursor.getSeqString() : null;
			if (sequence != null && !sequence.isEmpty()) {
				peptides.add(sequence);
			}
		}
		Map<String, Set<String>> peptideProteins = mapPeptideProteins(peptides, predictedFolder, allowedGenomes, 20);
		ArrayList<DiaNNPrecursor> result = new ArrayList<>();
		for (DiaNNPrecursor precursor : precursors) {
			if (precursor == null) {
				continue;
			}
			String sequence = precursor.getSeqString();
			if (sequence == null || sequence.isEmpty()) {
				continue;
			}
			Set<String> proteins = peptideProteins.get(sequence);
			if (proteins != null && !proteins.isEmpty()) {
				precursor.setProGroups(proteins.toArray(new String[0]));
				result.add(precursor);
			}
		}
		return result;
	}

	public static List<MetaPSM> mapPsmsToProteins(MetaPSM[] psms, File predictedFolder) {
		return mapPsmsToProteins(psms, predictedFolder, null);
	}

	public static List<MetaPSM> mapPsmsToProteins(MetaPSM[] psms, File predictedFolder, Set<String> allowedGenomes) {
		if (psms == null || psms.length == 0 || predictedFolder == null || !predictedFolder.isDirectory()) {
			return psms != null && psms.length > 0 ? new ArrayList<>(Arrays.asList(psms)) : new ArrayList<>();
		}
		HashSet<String> peptides = new HashSet<>();
		for (MetaPSM psm : psms) {
			String sequence = psm != null ? psm.getSequence() : null;
			if (sequence != null && !sequence.isEmpty()) {
				peptides.add(sequence);
			}
		}
		Map<String, Set<String>> peptideProteins = mapPeptideProteins(peptides, predictedFolder, allowedGenomes, 20);
		ArrayList<MetaPSM> result = new ArrayList<>();
		for (MetaPSM psm : psms) {
			if (psm == null) {
				continue;
			}
			String sequence = psm.getSequence();
			if (sequence == null || sequence.isEmpty()) {
				continue;
			}
			Set<String> proteins = peptideProteins.get(sequence);
			if (proteins != null && !proteins.isEmpty()) {
				psm.setProteins(proteins.toArray(new String[0]));
				result.add(psm);
			}
		}
		return result;
	}

	private static boolean isAllowedGenome(File file, Set<String> allowedGenomes) {
		if (allowedGenomes == null || allowedGenomes.isEmpty()) {
			return true;
		}
		String name = file.getName();
		String genome = name.substring(0, name.length() - PREDICTED_TSV_SUFFIX.length());
		return allowedGenomes.contains(genome);
	}

	private static void mapFile(File file, Set<String> peptides, Set<String> allowedProteins, Set<String> matchedPeptides,
			ConcurrentHashMap<String, HashMap<String, HashSet<String>>> genomeProteinPeptides) {
		try (BufferedReader reader = new BufferedReader(new FileReader(file))) {
			HashMap<String, HashSet<String>> proteinPeptides = new HashMap<>();
			reader.readLine();
			String line;
			while ((line = reader.readLine()) != null) {
				String[] columns = line.split("\t");
				if (columns.length < 2) {
					continue;
				}
				String protein = columns[0];
				String peptide = columns[1];
				if (!peptides.contains(peptide)) {
					continue;
				}
				if (allowedProteins != null && !allowedProteins.isEmpty() && !allowedProteins.contains(protein)) {
					continue;
				}
				matchedPeptides.add(peptide);
				proteinPeptides.computeIfAbsent(protein, ignored -> new HashSet<>()).add(peptide);
			}
			if (!proteinPeptides.isEmpty()) {
				String name = file.getName();
				String genome = name.substring(0, name.length() - PREDICTED_TSV_SUFFIX.length());
				genomeProteinPeptides.put(genome, proteinPeptides);
			}
		} catch (IOException e) {
			LOGGER.error("Error reading predicted peptide file {}", file.getAbsolutePath(), e);
		}
	}

	private static void mapProteinAssignments(File file, ConcurrentMap<String, Set<String>> peptideProteins) {
		try (BufferedReader reader = new BufferedReader(new FileReader(file))) {
			reader.readLine();
			String line;
			while ((line = reader.readLine()) != null) {
				String[] columns = line.split("\t");
				if (columns.length > 1 && peptideProteins.containsKey(columns[1])) {
					peptideProteins.get(columns[1]).add(columns[0]);
				}
			}
		} catch (IOException e) {
			LOGGER.warn("Error reading predicted TSV: {}", file.getAbsolutePath(), e);
		}
	}
}
