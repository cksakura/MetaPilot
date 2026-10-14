package uk.ac.quadram.metapilot.task.io.qpx;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Locale;

/**
 * One entry of the QPX {@code modifications} array: a modification (named when it
 * maps to UNIMOD, otherwise identified by its signed mass) and the positions it
 * occupies on the peptide.
 * <p>
 * Positions follow the QPX convention: 0 = peptide N-terminus, 1..n = residue
 * (1-based), n+1 = peptide C-terminus. A modification the engine reports without
 * a site (the unlocalized open-search delta mass) has {@code localized == false};
 * its positions list then holds the candidate sites, if any are known, and the
 * ProForma writer places it in the unlocalized {@code [mass]?} prefix.
 *
 * @see <a href="https://qpx.quantms.org/spec/modifications/">QPX modifications</a>
 */
public final class QpxModification {

	/** One site of a modification. {@code aminoAcid} is null for the termini. */
	public static final class Site {
		public final int position;
		public final String aminoAcid;

		public Site(int position, String aminoAcid) {
			this.position = position;
			this.aminoAcid = aminoAcid;
		}
	}

	private final String name;
	private final String accession;
	private final double mass;
	private final boolean localized;
	private final List<Site> sites;

	private QpxModification(String name, String accession, double mass, boolean localized, List<Site> sites) {
		this.name = name;
		this.accession = accession;
		this.mass = mass;
		this.localized = localized;
		this.sites = Collections.unmodifiableList(new ArrayList<>(sites));
	}

	/**
	 * A modification on one site of {@code sequence}. {@code position} uses the QPX
	 * numbering (0 = N-term, length+1 = C-term).
	 */
	public static QpxModification at(UnimodCatalog.Entry unimod, double massIfUnknown, String sequence,
			int position) {
		return at(unimod, massIfUnknown, null, sequence, position);
	}

	/**
	 * As {@link #at(UnimodCatalog.Entry, double, String, int)}, with the engine's own
	 * name as the last resort when neither a UNIMOD entry nor a mass is known.
	 */
	public static QpxModification at(UnimodCatalog.Entry unimod, double massIfUnknown, String engineName,
			String sequence, int position) {
		List<Site> sites = new ArrayList<>(1);
		sites.add(siteOf(sequence, position));
		return of(unimod, massIfUnknown, engineName, true, sites);
	}

	/** A modification whose site the engine did not determine. */
	public static QpxModification unlocalized(UnimodCatalog.Entry unimod, double massIfUnknown, String sequence,
			int[] candidatePositions) {
		List<Site> sites = new ArrayList<>();
		if (candidatePositions != null) {
			for (int p : candidatePositions) {
				sites.add(siteOf(sequence, p));
			}
		}
		return of(unimod, massIfUnknown, null, false, sites);
	}

	private static QpxModification of(UnimodCatalog.Entry unimod, double mass, String engineName,
			boolean localized, List<Site> sites) {
		if (unimod != null) {
			return new QpxModification(unimod.title, unimod.accession(), unimod.monoMass, localized, sites);
		}
		// QPX: "When no UNIMOD mapping exists the name field carries the mass with sign
		// and the accession field is null." Only when the engine gave a name but no
		// mass (and the name is not a UNIMOD title) is the engine's name kept instead.
		if (Double.isNaN(mass) && engineName != null && !engineName.isEmpty()) {
			return new QpxModification(engineName, null, mass, localized, sites);
		}
		return new QpxModification(formatMass(mass), null, mass, localized, sites);
	}

	static Site siteOf(String sequence, int position) {
		int length = sequence == null ? 0 : sequence.length();
		if (position <= 0 || position > length) {
			return new Site(position <= 0 ? 0 : length + 1, null);
		}
		return new Site(position, String.valueOf(sequence.charAt(position - 1)));
	}

	/** Signed mass with four decimals, e.g. {@code +79.9663}, the QPX/ProForma form. */
	public static String formatMass(double mass) {
		return String.format(Locale.ROOT, "%+.4f", mass);
	}

	public String getName() {
		return name;
	}

	/** {@code UNIMOD:n}, or null when the modification has no UNIMOD entry. */
	public String getAccession() {
		return accession;
	}

	/** Monoisotopic mass shift in Da (UNIMOD's when mapped, the engine's otherwise). */
	public double getMass() {
		return mass;
	}

	public boolean isLocalized() {
		return localized;
	}

	public List<Site> getSites() {
		return sites;
	}
}
