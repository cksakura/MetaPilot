package uk.ac.quadram.metapilot.core.msbooster;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

/**
 * Reads and subsets Percolator {@code .pin} files, the format MSBooster consumes
 * and enriches.
 *
 * <h2>The pin is not narrowed before prediction</h2>
 * MSBooster gets the pin exactly as MSFragger wrote it, which is what FragPipe
 * does. Narrowing it to the PSMs stage 1 consumes looks like a saving, but a
 * narrowed pin cannot keep the {@code <sample>.pin} base name MSBooster pairs
 * spectra on, and MSBooster already fits its retention-time calibration on a
 * confident subset it selects itself. The narrowing happens here instead, on the
 * way back out: {@link #readFeatures} keys every row, and a caller looks up only
 * the PSMs it cares about, so unused rows cost nothing beyond the read.
 *
 * <h2>Keying</h2>
 * A pin row is identified by MSFragger's {@code SpecId}, shaped
 * {@code <basename>.<scan>.<scan>.<charge>_<rank>}. That is decomposed into
 * {@code basename:scan:rank}, which is what an {@code OpenSearchPsm} carries as
 * {@code getFileName()}, {@code getScan()} and {@code getHitRank()}. Matching on
 * the parsed triple rather than the raw string means a pin whose SpecId
 * formatting differs slightly - zero padding, a missing rank suffix - still
 * lines up instead of silently matching nothing.
 */
public final class MSBoosterPin {

	private static final String SPEC_ID = "SpecId";
	private static final String DEFAULT_DIRECTION = "DefaultDirection";

	private MSBoosterPin() {
	}

	/**
	 * Canonical key for a PSM: {@code file:scan:rank}. Scan is parsed to an int so
	 * zero padding does not matter.
	 */
	public static String key(String fileBaseName, int scan, int hitRank) {
		return fileBaseName + ":" + scan + ":" + hitRank;
	}

	/**
	 * Parses a pin {@code SpecId} into the same key. Returns null when the shape
	 * is not recognised, so a caller can count and report unmatched rows rather
	 * than mis-associating features.
	 */
	public static String keyFromSpecId(String specId) {
		if (specId == null || specId.isBlank()) {
			return null;
		}
		String id = specId.trim();
		int rank = 1;
		int underscore = id.lastIndexOf('_');
		if (underscore > 0) {
			String tail = id.substring(underscore + 1);
			try {
				rank = Integer.parseInt(tail);
				id = id.substring(0, underscore);
			} catch (NumberFormatException ignored) {
				// No rank suffix; the whole string is the spectrum id.
			}
		}
		String[] parts = id.split("\\.");
		if (parts.length < 4) {
			return null;
		}
		// Last three are scan, scan, charge; everything before them is the base
		// name, which may itself contain dots.
		StringBuilder base = new StringBuilder();
		for (int i = 0; i < parts.length - 3; i++) {
			if (i > 0) {
				base.append('.');
			}
			base.append(parts[i]);
		}
		try {
			int scan = Integer.parseInt(parts[parts.length - 3]);
			return key(base.toString(), scan, rank);
		} catch (NumberFormatException e) {
			return null;
		}
	}

	/** What {@link #readFeatures} found. */
	public static final class Features {
		/** Feature columns actually present, in pin order. */
		public final List<String> names;
		/** PSM key to that row's values, aligned with {@link #names}. */
		public final Map<String, double[]> byKey;
		/** Requested columns the file did not contain; empty when all were found. */
		public final List<String> missing;

		Features(List<String> names, Map<String, double[]> byKey, List<String> missing) {
			this.names = names;
			this.byKey = byKey;
			this.missing = missing;
		}

		public boolean isEmpty() {
			return names.isEmpty() || byKey.isEmpty();
		}
	}

	/**
	 * Reads the named columns out of an edited pin.
	 * <p>
	 * Requesting columns by name rather than diffing against the input pin means
	 * a column MSBooster renames upstream is reported as missing and logged,
	 * instead of quietly changing what the classifier trains on.
	 */
	public static Features readFeatures(File editedPin, Collection<String> wanted) throws IOException {
		if (editedPin == null || !editedPin.isFile()) {
			throw new IOException("Edited .pin not found: " + editedPin);
		}
		Map<String, double[]> byKey = new HashMap<>();
		List<String> present = new ArrayList<>();
		List<String> missing = new ArrayList<>();
		Set<String> requested = new LinkedHashSet<>(wanted);

		try (BufferedReader in = Files.newBufferedReader(editedPin.toPath(), StandardCharsets.UTF_8)) {
			String header = in.readLine();
			if (header == null) {
				throw new IOException("Empty edited .pin: " + editedPin);
			}
			String[] cols = header.split("\t", -1);
			Map<String, Integer> index = new LinkedHashMap<>();
			for (int i = 0; i < cols.length; i++) {
				index.put(cols[i].trim(), i);
			}
			if (!index.containsKey(SPEC_ID) && !cols[0].trim().equalsIgnoreCase(SPEC_ID)) {
				throw new IOException("Not a .pin file (no SpecId column): " + editedPin);
			}
			List<Integer> take = new ArrayList<>();
			for (String name : requested) {
				Integer at = index.get(name);
				if (at == null) {
					missing.add(name);
				} else {
					present.add(name);
					take.add(at);
				}
			}
			if (present.isEmpty()) {
				return new Features(present, byKey, missing);
			}
			String line;
			boolean first = true;
			while ((line = in.readLine()) != null) {
				if (line.isBlank()) {
					continue;
				}
				if (first && line.startsWith(DEFAULT_DIRECTION)) {
					first = false;
					continue;
				}
				first = false;
				String[] f = line.split("\t", -1);
				String k = keyFromSpecId(f[0]);
				if (k == null) {
					continue;
				}
				double[] values = new double[take.size()];
				boolean usable = true;
				for (int i = 0; i < take.size(); i++) {
					int at = take.get(i);
					if (at >= f.length) {
						usable = false;
						break;
					}
					try {
						values[i] = Double.parseDouble(f[at].trim());
					} catch (NumberFormatException e) {
						// A non-numeric cell is left as a missing value rather than
						// zero: zero is a meaningful similarity, absent is not.
						values[i] = Double.NaN;
					}
				}
				if (usable) {
					byKey.put(k, values);
				}
			}
		}
		if (!missing.isEmpty()) {
			// Not a warning. A classifier silently trained on five of seven
			// predicted features produces a plausible-looking result that cannot
			// be compared with any other run, and nothing downstream records which
			// features were present. If MSBooster renamed or dropped a column, that
			// is a version mismatch to be fixed, not absorbed.
			throw new IOException("Edited pin " + editedPin.getName() + " is missing MSBooster feature column(s) "
					+ missing + "; it has " + present + ". MSBooster 1.4.x writes all of "
					+ MSBoosterHeadlessRunner.DEFAULT_FEATURE_COLUMNS + " - a different set means a different "
					+ "MSBooster version or a half-written file.");
		}
		return new Features(present, byKey, missing);
	}
}
