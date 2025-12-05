/**
 * 
 */
package uk.ac.quadram.metapilot.dbSearch.philosopher;

import java.io.File;

/**
 * @author Kai Cheng
 *
 */
public class PhilosopherPepReader {
	
	public PhilosopherPepReader(String in) {
		this(new File(in));
	}
	
	public PhilosopherPepReader(File in) {
		this.read(in);
	}
	
	private void read(File in) {

	}
	
}
