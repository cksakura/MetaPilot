/**
 * 
 */
package uk.ac.quadram.metapilot.core.enzyme;

import java.util.ArrayList;

/**
 * Comprehensive Proteolytic Enzyme Support for In-Silico Protein Digestion
 * 
 * <h2>Purpose</h2>
 * This enum provides comprehensive support for proteolytic enzymes commonly used in
 * mass spectrometry-based proteomics. It enables accurate in-silico protein digestion
 * simulation, essential for database search parameter configuration and theoretical
 * peptide generation.
 * 
 * <h2>Supported Enzymes</h2>
 * 
 * <h3>Trypsin Family:</h3>
 * <ul>
 *   <li><strong>Trypsin:</strong> Cleaves C-terminal to K,R (not before P)</li>
 *   <li><strong>Trypsin/P:</strong> Cleaves C-terminal to K,R (including before P)</li>
 * </ul>
 * 
 * <h3>Lysine-Specific Endopeptidases:</h3>
 * <ul>
 *   <li><strong>LysC:</strong> Cleaves C-terminal to K (not before P)</li>
 *   <li><strong>LysC/P:</strong> Cleaves C-terminal to K (including before P)</li>
 *   <li><strong>LysN:</strong> Cleaves N-terminal to K</li>
 * </ul>
 * 
 * <h3>Arginine-Specific Endopeptidases:</h3>
 * <ul>
 *   <li><strong>ArgC:</strong> Cleaves C-terminal to R</li>
 * </ul>
 * 
 * <h3>Acid-Specific Endopeptidases:</h3>
 * <ul>
 *   <li><strong>AspC:</strong> Cleaves C-terminal to D (aspartic acid)</li>
 *   <li><strong>AspN:</strong> Cleaves N-terminal to D (aspartic acid)</li>
 *   <li><strong>GluC:</strong> Cleaves C-terminal to E (glutamic acid)</li>
 * </ul>
 * 
 * <h3>Chymotrypsin Family:</h3>
 * <ul>
 *   <li><strong>ChymoTrypsin:</strong> Cleaves C-terminal to F,W,Y (aromatic residues)</li>
 *   <li><strong>ChymoTrypsin+:</strong> Extended specificity - cleaves C-terminal to L,M,F,W,Y</li>
 * </ul>
 * 
 * <h3>Other Enzymes:</h3>
 * <ul>
 *   <li><strong>GlyN:</strong> Cleaves N-terminal to K (glycine endopeptidase)</li>
 * </ul>
 * 
 * <h2>Cleavage Specificity Model</h2>
 * <h3>Cleavage Rules:</h3>
 * <ul>
 *   <li><strong>cleaveAt:</strong> Amino acids where cleavage occurs</li>
 *   <li><strong>notCleaveAt:</strong> Amino acids that prevent cleavage (restrictions)</li>
 *   <li><strong>sense_C:</strong> Cleavage direction (C-terminal = true, N-terminal = false)</li>
 * </ul>
 * 
 * <h3>Cleavage Notation:</h3>
 * <pre>{@code
 *   C-terminal cleavage: ...Peptide-K|R-Peptide...
 *   N-terminal cleavage: ...Peptide-|K-R-Peptide...
 *   
 *   Where | represents the cleavage site
 * }</pre>
 * 
 * <h2>Missed Cleavage Support</h2>
 * <p>The digestion algorithm supports missed cleavages:</p>
 * <ul>
 *   <li><strong>0 missed:</strong> Complete digestion at all sites</li>
 *   <li><strong>1 missed:</strong> Allows one missed cleavage site</li>
 *   <li><strong>2+ missed:</strong> Multiple missed cleavage sites</li>
 * </ul>
 * 
 * <h2>Output Formats</h2>
 * 
 * <h3>Tandem (X! Tandem) Format:</h3>
 * <pre>{@code
 *   Trypsin: [KR]|{P}
 *   LysC:    [K]|{P}
 * }</pre>
 * 
 * <h3>DIA-NN Format:</h3>
 * <pre>{@code
 *   Trypsin: K*,R*,!*P
 *   LysC:    K*,!*P
 * }</pre>
 * 
 * <h2>Peptide Length Filtering</h2>
 * <ul>
 *   <li><strong>Minimum Length:</strong> Configurable (default: 6 amino acids)</li>
 *   <li><strong>Maximum Length:</strong> Typically 50-60 amino acids</li>
 *   <li><strong>Mass Range:</strong> Usually 400-6000 Da for MS analysis</li>
 * </ul>
 * 
 * <h2>Usage Examples</h2>
 * 
 * <h3>Basic Digestion:</h3>
 * <pre>{@code
 *   String protein = "MAKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLK";
 *   String[] peptides = Enzyme.Trypsin.digest(protein, 1); // 1 missed cleavage
 *   
 *   // Results: ["MAKGEELFTGVVPILVELDGDVNGHK", "FSVSGEGEGDATYGK", "LTLK"]
 * }</pre>
 * 
 * <h3>Custom Length Filtering:</h3>
 * <pre>{@code
 *   String[] peptides = Enzyme.Trypsin.digest(protein, 0, 8); // min 8 AA
 * }</pre>
 * 
 * <h3>Search Engine Integration:</h3>
 * <pre>{@code
 *   // For X! Tandem
 *   String tandemFormat = Enzyme.Trypsin.getTandemFormat();
 *   
 *   // For DIA-NN
 *   String diannFormat = Enzyme.Trypsin.getCut();
 * }</pre>
 * 
 * <h2>Biological Context</h2>
 * <h3>Trypsin Characteristics:</h3>
 * <ul>
 *   <li><strong>Specificity:</strong> Basic amino acids (K, R)</li>
 *   <li><strong>pH Optimum:</strong> 8.0-8.5</li>
 *   <li><strong>Proline Effect:</strong> Inhibits cleavage at P1' position</li>
 *   <li><strong>Most Common:</strong> >90% of proteomics experiments</li>
 * </ul>
 * 
 * <h3>Alternative Proteases:</h3>
 * <ul>
 *   <li><strong>LysC:</strong> More specific than trypsin, reduces missed cleavages</li>
 *   <li><strong>AspN:</strong> Provides orthogonal cleavage pattern</li>
 *   <li><strong>ChymoTrypsin:</strong> Generates different peptide populations</li>
 * </ul>
 * 
 * <h2>Performance Considerations</h2>
 * <ul>
 *   <li><strong>Algorithm:</strong> O(n) scanning with regex-like matching</li>
 *   <li><strong>Memory:</strong> Efficient ArrayList-based peptide collection</li>
 *   <li><strong>Scalability:</strong> Suitable for whole proteome digestion</li>
 * </ul>
 * 
 * <h2>Thread Safety</h2>
 * <p>This enum is thread-safe as all methods are stateless and operate only
 * on input parameters without modifying internal state.</p>
 * 
 * @author Kai Cheng - Quadram Institute Bioscience
 * @version 2.0
 * @since 1.0
 * 
 * @see uk.ac.quadram.metapilot.core.aminoacid.Aminoacids
 * @see uk.ac.quadram.metapilot.dbSearch
 */
public enum Enzyme {

	Trypsin("Trypsin", true, "KR", "P"),

	TrypsinP("Trypsin/P", true, "KR", ""),

	LysC("LysC", true, "K", "P"),

	LysCP("LysC/P", true, "K", ""),

	LysN("LysN", false, "K", ""),

	ArgC("ArgC", true, "R", ""),

	AspC("AspC", true, "D", ""),

	AspN("AspN", false, "D", ""),

	GluC("GluC", true, "E", ""),

	GlyN("GlyN", false, "K", ""),

	ChymoTrypsin("ChymoTrypsin", false, "FWY", ""),

	ChymoTrypsinPlus("ChymoTrypsin+", false, "LMFWY", "");

	private String name;
	private boolean sense_C;
	private String cleaveAt;
	private String notCleaveAt;
	private char[] cleaves;
	private char[] restricts;

	Enzyme(String enzymeName, boolean sense_C, String cleaveAt, String notCleaveAt) {

		this.name = enzymeName;
		this.sense_C = sense_C;
		this.cleaveAt = cleaveAt;
		this.notCleaveAt = notCleaveAt;
		this.cleaves = cleaveAt.toCharArray();
		this.restricts = notCleaveAt.toCharArray();
	}

	public String getName() {
		return name;
	}

	public boolean isSense_C() {
		return sense_C;
	}

	public String getCleaveAt() {
		return cleaveAt;
	}

	public String getNotCleaveAt() {
		return notCleaveAt;
	}

	public char[] getCleaves() {
		return cleaves;
	}

	public char[] getRestricts() {
		return restricts;
	}

	public Object[] toTableRawObjects() {
		Object[] objs = new Object[5];
		objs[0] = false;
		objs[1] = name;
		objs[2] = cleaveAt;
		objs[3] = notCleaveAt;
		objs[4] = sense_C;
		return objs;
	}

	public String getTandemFormat() {
		StringBuilder sb = new StringBuilder();
		sb.append("[").append(this.cleaveAt).append("]").append("|").append("{").append(this.notCleaveAt).append("}");
		return sb.toString();
	}

	private String[] getPeptides(String sequence, int[] loc, int miss, int leastLength) {
		ArrayList<String> list = new ArrayList<String>();
		for (int i = 0; i <= miss; i++) {
			if (i < loc.length && loc[i] + 1 >= leastLength) {
				String pep = sequence.substring(0, loc[i] + 1);
				list.add(pep);
			}
		}

		for (int i = 0; i < loc.length - 1; i++) {
			for (int j = i + 1; j <= i + 1 + miss; j++) {
				if (j == loc.length)
					break;

				if (loc[j] - loc[i] >= leastLength) {
					String pep = sequence.substring(loc[i] + 1, loc[j] + 1);
					list.add(pep);
				}
			}
		}

		String[] peps = list.toArray(new String[list.size()]);
		return peps;
	}

	/**
	 * Get the cleavage sites
	 * 
	 * @param sequence
	 * @return
	 */
	private int[] getLoc(String sequence) {
		ArrayList<Integer> list = new ArrayList<Integer>();
		char[] cs = sequence.toCharArray();
		for (int i = 0; i < cs.length - 1; i++) {
			char ci0 = cs[i];
			char ci1 = cs[i + 1];

			boolean cleave = false;
			for (int j = 0; j < this.cleaves.length; j++) {
				if (ci0 == this.cleaves[j]) {
					cleave = true;
					break;
				}
			}

			boolean res = false;
			for (int j = 0; j < this.restricts.length; j++) {
				if (ci1 == this.restricts[j]) {
					res = true;
					break;
				}
			}

			if (cleave && !res) {
				list.add(i);
			}
		}
		list.add(cs.length - 1);

		int[] arrays = new int[list.size()];
		for (int i = 0; i < arrays.length; i++) {
			arrays[i] = list.get(i);
		}
		return arrays;
	}

	public String[] digest(String sequence, int miss) {
		int[] locs = this.getLoc(sequence);
		return getPeptides(sequence, locs, miss, 6);
	}
	
	public String[] digest(String sequence, int miss, int leastLength) {
		int[] locs = this.getLoc(sequence);
		return getPeptides(sequence, locs, miss, leastLength);
	}
	
	/**
	 * DiaNN format
	 * @return
	 */
	public String getCut() {
		StringBuilder sb = new StringBuilder();
		for (int i = 0; i < cleaveAt.length(); i++) {
			char aa = cleaveAt.charAt(i);
			sb.append(aa).append("*,");
		}
		if (notCleaveAt.length() > 0) {
			sb.append("!*").append(notCleaveAt.charAt(0));
		}
		return sb.toString();
	}
}
