package uk.ac.quadram.metapilot.task.io.qpx;

import java.util.ArrayList;
import java.util.List;
import java.util.Locale;

/**
 * Builds ProForma 2.0 peptidoform strings from a base sequence and QPX
 * modifications, following the QPX conventions: a UNIMOD accession when the
 * modification has one ({@code PEPT[UNIMOD:21]IDE}), otherwise the signed mass
 * ({@code PEPT[+79.9663]IDE}); N-terminal {@code [..]-PEP}, C-terminal
 * {@code PEP-[..]}; and a modification without a site in the unlocalized prefix
 * ({@code [+14.0157]?PEPTIDE}).
 *
 * @see <a href="https://qpx.quantms.org/spec/modifications/">QPX modifications</a>
 */
public final class ProForma {

	private ProForma() {
	}

	public static String build(String sequence, List<QpxModification> modifications) {
		String seq = sequence == null ? "" : sequence.toUpperCase(Locale.ROOT);
		if (modifications == null || modifications.isEmpty()) {
			return seq;
		}
		int n = seq.length();
		List<List<String>> residueTags = new ArrayList<>(n + 2);
		for (int i = 0; i < n + 2; i++) {
			residueTags.add(null);
		}
		StringBuilder unlocalized = new StringBuilder();
		for (QpxModification mod : modifications) {
			String tag = "[" + tokenOf(mod) + "]";
			if (!mod.isLocalized() || mod.getSites().isEmpty()) {
				unlocalized.append(tag).append('?');
				continue;
			}
			for (QpxModification.Site site : mod.getSites()) {
				int p = Math.max(0, Math.min(site.position, n + 1));
				List<String> tags = residueTags.get(p);
				if (tags == null) {
					tags = new ArrayList<>(1);
					residueTags.set(p, tags);
				}
				tags.add(tag);
			}
		}

		StringBuilder sb = new StringBuilder(n + 32);
		sb.append(unlocalized);
		List<String> nTerm = residueTags.get(0);
		if (nTerm != null) {
			for (String t : nTerm) {
				sb.append(t);
			}
			sb.append('-');
		}
		for (int i = 1; i <= n; i++) {
			sb.append(seq.charAt(i - 1));
			List<String> tags = residueTags.get(i);
			if (tags != null) {
				for (String t : tags) {
					sb.append(t);
				}
			}
		}
		List<String> cTerm = residueTags.get(n + 1);
		if (cTerm != null) {
			sb.append('-');
			for (String t : cTerm) {
				sb.append(t);
			}
		}
		return sb.toString();
	}

	/**
	 * Accession when mapped, else the name - which for an unmapped modification is
	 * already its signed mass (see {@link QpxModification}).
	 */
	static String tokenOf(QpxModification mod) {
		if (mod.getAccession() != null) {
			return mod.getAccession();
		}
		return mod.getName();
	}
}
