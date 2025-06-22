package uk.ac.quadram.metapilot.task;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.PrintWriter;
import java.util.Date;

import org.apache.commons.io.FileUtils;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import uk.ac.quadram.metapilot.gui.JavaFXProgressBar;
import uk.ac.quadram.metapilot.resources.MetaResourceParameters;
import uk.ac.quadram.metapilot.task.par.MetaBaseParameter;

/**
 * Create the sample-specific protein database. <br>
 * Two workflow "spectra cluster" and "MetaProIQ" are provided in this step.
 * 
 * @author Kai Cheng
 *
 */
public abstract class MetaDbCreateTask extends MetaAbstractTask {

	protected File dbcFile;
	protected File ssdbFile;
	protected String originalFasta;
	protected File[] clusterResults;

	protected Logger LOGGER = LogManager.getLogger(MetaDbCreateTask.class);
	protected String taskName = "Sample specific database generation";

	public MetaDbCreateTask(MetaBaseParameter metaPar, MetaResourceParameters mrp, JavaFXProgressBar progressBar) {
		super(metaPar, mrp, progressBar);
	}

	@Override
	protected void initial() {
		// TODO Auto-generated method stub

	}

	@Override
	protected Boolean call() throws Exception {
		// Update progress at start
		updateProgress(0, 100);
		updateMessage("Starting database creation...");

		progressBar.setString(taskName);

		LOGGER.info(taskName + ": start");
		System.out.println(format.format(new Date()) + "\t" + taskName + ": start");

		if (metaPar.isSpectraCluster()) {
			this.dbcFile = new File(metaPar.getResult(), "db_generation_msCluster");
		} else {
			this.dbcFile = new File(metaPar.getResult(), "db_generation_metaProIq");
		}
		if (!dbcFile.exists()) {
			dbcFile.mkdir();
		}

		this.ssdbFile = new File(dbcFile, "sample_specific_db.fasta");

		if (this.ssdbFile.exists()) {

			if (ssdbFile.length() > 0) {
				LOGGER.info(taskName + ": sample specific database was found in " + ssdbFile + ", go to next step");
				System.out.println(format.format(new Date()) + "\t" + taskName
						+ ": sample specific database was found in " + ssdbFile + ", go to next step");
				metaPar.setCurrentDb(ssdbFile.getAbsolutePath());

				setProgress(35);

				return true;

			} else {

				FileUtils.deleteQuietly(ssdbFile);
			}
		}

		if (metaPar.isAppendHostDb()) {

			File combinedFile = new File(dbcFile, "combined_db.fasta");
			if (combinedFile.exists() && combinedFile.length() == 0) {

				FileUtils.deleteQuietly(combinedFile);
			}

			if (!combinedFile.exists()) {
				String microDb = metaPar.getMicroDb();
				String hostDb = metaPar.getHostDB();

				PrintWriter writer = new PrintWriter(combinedFile);

				BufferedReader r0 = new BufferedReader(new FileReader(microDb));
				String l0 = null;
				while ((l0 = r0.readLine()) != null) {
					writer.println(l0);
				}
				r0.close();

				BufferedReader r1 = new BufferedReader(new FileReader(hostDb));
				String l1 = null;
				while ((l1 = r1.readLine()) != null) {
					writer.println(l1);
				}
				r1.close();

				writer.close();
			}

			metaPar.setCurrentDb(combinedFile.getAbsolutePath());
			originalFasta = combinedFile.getAbsolutePath();

		} else {
			metaPar.setCurrentDb(metaPar.getMicroDb());
			originalFasta = metaPar.getMicroDb();
		}

		int proCount = 0;

		if (metaPar.isSpectraCluster()) {

			cluster();

			clusterResults = clusterSearch();

			proCount = writeDatabase(true);

		} else {

			clusterResults = metaProIQ();

			proCount = writeDatabase(false);
		}

		LOGGER.info(taskName + ": sample specific protein database " + ssdbFile + " was created, " + proCount
				+ " target protein sequences were included");
		System.out.println(format.format(new Date()) + "\t" + taskName + ": sample specific protein database "
				+ ssdbFile + " was created, " + proCount + " target protein sequences were included");

		metaPar.setCurrentDb(ssdbFile.getAbsolutePath());

		if (proCount == 0) {
			System.err.println(format.format(new Date()) + "\t" + taskName
					+ ": cannot create the sample-specific database because no proteins are identified from the original database, task terminated.");
			LOGGER.error(taskName
					+ ": cannot create the sample-specific database because no proteins are identified from the original database, task terminated.");
			return false;
		}

		setProgress(35);

		LOGGER.info(taskName + ": finished");
		System.out.println(format.format(new Date()) + "\t" + taskName + ": finished");

		return ssdbFile.exists();
	}

	/**
	 * Clustering strategy
	 * 
	 */
	protected abstract void cluster();

	/**
	 * Clustering strategy
	 * 
	 * 
	 */
	protected abstract File[] clusterSearch();

	/**
	 * MetaProIQ strategy
	 */
	protected abstract File[] metaProIQ();

	/**
	 * Write the protein sequences into the sample-specific database
	 * 
	 * @return
	 */
	protected abstract int writeDatabase(boolean cluster);

	@Override
	protected String getTaskName() {
		// TODO Auto-generated method stub
		return taskName;
	}

	@Override
	protected Logger getLogger() {
		// TODO Auto-generated method stub
		return LOGGER;
	}

	public void forceStop() {
		this.cancel();
	}
}
