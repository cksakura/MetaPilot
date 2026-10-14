/**
 * 
 */
package uk.ac.quadram.metapilot.dbSearch.pfind;

import java.util.HashMap;

import uk.ac.quadram.metapilot.core.mod.PostTransModification;
import uk.ac.quadram.metapilot.task.io.pep.MetaPSM;

/**
 * @author Kai Cheng
 *
 */
public class PFindPSM extends MetaPSM {

	private double expMH;
	private double calMH;
	private double massShift;
	private double rawScore;
	private double finalScore;
	private String mod;
	private PostTransModification ptm;
	private int specificity;
	private boolean isTarget;
	private int miss;
	private double fragMassShift;
	private double rt;

	public PFindPSM(String fileName, int scan, double expMH, int charge, double qValue, String sequence, double calMH,
			double massShift, double rawScore, double finalScore, String mod, PostTransModification ptm,
			int specificity, String[] proteins, boolean isTarget, int miss, double fragMassShift) {
		super(fileName, scan, sequence, buildModSeq(sequence, mod), charge, proteins, finalScore, qValue);
		this.expMH = expMH;
		this.calMH = calMH;
		this.massShift = massShift;
		this.rawScore = rawScore;
		this.finalScore = finalScore;
		this.mod = mod;
		this.ptm = ptm;
		this.specificity = specificity;
		this.isTarget = isTarget;
		this.miss = miss;
		this.fragMassShift = fragMassShift;
	}

	private static final java.util.regex.Pattern PFIND_SITE = java.util.regex.Pattern
			.compile("[A-Z]|AnyN-term|AnyC-term|ProteinN-term|ProteinC-term");

	/**
	 * The modification name of one pFind "Modification" entry without its site:
	 * {@code Oxidation[M]} -> {@code Oxidation}, {@code Cation_Fe[II][D]} ->
	 * {@code Cation_Fe[II]}, and for an entry carrying pFind's alias group,
	 * {@code Methyl[S](Ser->Thr[S])} -> {@code Methyl}. The alias case used to be cut
	 * at the last '[' - inside the alias - giving the unbalanced
	 * {@code Methyl[S](Ser->Thr} (1,501 of 1,255,000 PSMs in the 12mix pFind psms.tsv).
	 */
	static String modificationName(String entry) {
		String e = entry.trim();
		if (e.endsWith(")")) {
			int k = e.indexOf("](");
			while (k >= 0) {
				String primary = e.substring(0, k + 1);
				int open = primary.lastIndexOf('[');
				if (open > 0 && PFIND_SITE.matcher(primary.substring(open + 1, primary.length() - 1)).matches()) {
					e = primary;
					break;
				}
				k = e.indexOf("](", k + 1);
			}
		}
		int bracket = e.lastIndexOf('[');
		return bracket > 0 ? e.substring(0, bracket) : e;
	}

	private static String buildModSeq(String sequence, String mod) {
		if (sequence == null) {
			return "";
		}
		if (mod == null || mod.trim().isEmpty()) {
			return sequence;
		}

		HashMap<Integer, String> modmap = new HashMap<>();
		StringBuilder sb = new StringBuilder();
		String[] mods = mod.split(";");
		for (int i = 0; i < mods.length; i++) {
			int id = mods[i].indexOf(",");
			if (id > 0) {
				int loc = Integer.parseInt(mods[i].substring(0, id));
				String name = modificationName(mods[i].substring(id + 1));
				modmap.put(loc, name);
			}
		}

		for (int i = 0; i < sequence.length(); i++) {
			if (modmap.containsKey(i)) {
				sb.append("(").append(modmap.get(i)).append(")");
			}
			sb.append(sequence.charAt(i));
		}
		if (modmap.containsKey(sequence.length())) {
			sb.append("(").append(modmap.get(sequence.length())).append(")");
		}
		if (modmap.containsKey(sequence.length() + 1)) {
			sb.append("(").append(modmap.get(sequence.length() + 1)).append(")");
		}
		return sb.toString();
	}

	public double getExpMH() {
		return expMH;
	}

	public double getCalMH() {
		return calMH;
	}

	public double getMassShift() {
		return massShift;
	}

	public double getRawScore() {
		return rawScore;
	}

	public double getFinalScore() {
		return finalScore;
	}

	public String getMod() {
		return mod;
	}

	public PostTransModification getPtm() {
		return ptm;
	}

	public int getSpecificity() {
		return specificity;
	}

	public void setProteins(String[] proteins) {
		this.proteins = proteins;
	}

	public String getProtein() {
		StringBuilder sb = new StringBuilder();
		for (String pro : proteins) {
			sb.append(pro).append(";");
		}
		if (sb.length() > 0) {
			return sb.substring(0, sb.length() - 1);
		} else {
			return "";
		}
	}

	public boolean isTarget() {
		return isTarget;
	}

	@Override
	public int getMiss() {
		return miss;
	}

	public double getFragMassShift() {
		return fragMassShift;
	}

	@Override
	public double getRt() {
		return rt;
	}

	public void setRt(double rt) {
		this.rt = rt;
	}

	public String getModseq() {
		return buildModSeq(baseSequence, mod);
	}

	/**
	 * No isobaric tags
	 * 
	 * @return
	 */
	public String getModseqNoIsoTags() {

		if (mod.trim().length() == 0) {
			return baseSequence;
		}

		HashMap<Integer, String> modmap = new HashMap<Integer, String>();
		StringBuilder sb = new StringBuilder();
		String[] mods = mod.split(";");
		for (int i = 0; i < mods.length; i++) {
			int id = mods[i].indexOf(",");
			if (id > 0) {
				int loc = Integer.parseInt(mods[i].substring(0, id));
				String name = mods[i].substring(id + 1, mods[i].indexOf("["));

				if (!name.contains("plex")) {
					modmap.put(loc, name);
				}
			}
		}

		for (int i = 0; i < baseSequence.length(); i++) {
			if (modmap.containsKey(i)) {
				sb.append("[").append(modmap.get(i)).append("]");
			}
			sb.append(baseSequence.charAt(i));
		}
		if (modmap.containsKey(baseSequence.length())) {
			sb.append("[").append(modmap.get(baseSequence.length())).append("]");
		}
		if (modmap.containsKey(baseSequence.length() + 1)) {
			sb.append("[").append(modmap.get(baseSequence.length() + 1)).append("]");
		}
		return sb.toString();
	}

	@Override
	public double getExpMass() {
		return expMH;
	}

	@Override
	public double getCalMass() {
		return calMH;
	}

	@Override
	public double getMassDiff() {
		return massShift;
	}

}
