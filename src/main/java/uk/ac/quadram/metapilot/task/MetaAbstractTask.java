package uk.ac.quadram.metapilot.task;

import java.text.SimpleDateFormat;
import java.util.Date;
import java.util.concurrent.ExecutionException;

import javafx.application.Platform;
import javafx.concurrent.Task;

import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.core.config.builder.api.ConfigurationBuilder;
import org.apache.logging.log4j.core.config.builder.impl.BuiltConfiguration;

import uk.ac.quadram.metapilot.gui.JavaFXProgressBar;
import uk.ac.quadram.metapilot.resources.MetaResourceParameters;
import uk.ac.quadram.metapilot.task.par.MetaBaseParameter;

/**
 * An abstract class for background tasks in MetaPilot that provides progress
 * tracking, lifecycle management and backward compatibility with the
 * SwingWorker-based design.
 * 
 * This class extends JavaFX's Task framework while maintaining an API
 * compatible with the original Swing implementation.
 * 
 * @author Kai Cheng
 */
public abstract class MetaAbstractTask extends Task<Boolean> {

	protected MetaBaseParameter metaPar;
	protected MetaResourceParameters mrp;
	protected JavaFXProgressBar progressBar;
	protected Task<Boolean> nextWork;
	protected SimpleDateFormat format = new SimpleDateFormat("MM_dd_yyyy_HH_mm_ss");
	protected ConfigurationBuilder<BuiltConfiguration> builder;

	/**
	 * Create a new MetaAbstractTask with the given parameters
	 * 
	 * @param metaPar      The meta parameter
	 * @param msv          The meta sources
	 * @param progressBar1 The primary progress bar
	 * @param progressBar2 The secondary progress bar
	 * @param nextWork     The next task to execute after this one completes
	 */
	public MetaAbstractTask(MetaBaseParameter metaPar, MetaResourceParameters mrp, JavaFXProgressBar progressBar) {
		this.metaPar = metaPar;
		this.mrp = mrp;
		this.progressBar = progressBar;
		// Set up progress property listeners
		this.progressProperty().addListener((observable, oldValue, newValue) -> {
			if (this.progressBar != null) {
				this.progressBar.setValue(newValue.intValue() * 100);
			}
		});

		// Set up handlers for task completion
		this.setOnSucceeded(event -> done());
		this.setOnFailed(event -> done());
		this.setOnCancelled(event -> {
			Logger LOGGER = getLogger();
			String taskName = getTaskName();
			LOGGER.info(taskName + ": cancelled");
			System.out.println(format.format(new Date()) + "\t" + taskName + ": task cancelled");
		});

		initial();
	}

	/**
	 * Create a new MetaAbstractTask with the given parameters
	 * 
	 * @param metaPar      The meta parameter
	 * @param msv          The meta sources
	 * @param progressBar1 The primary progress bar
	 * @param progressBar2 The secondary progress bar
	 * @param nextWork     The next task to execute after this one completes
	 */
	public MetaAbstractTask(MetaBaseParameter metaPar, MetaResourceParameters mrp, JavaFXProgressBar progressBar,
			Task<Boolean> nextWork) {
		this.metaPar = metaPar;
		this.mrp = mrp;
		this.progressBar = progressBar;
		this.nextWork = nextWork;

		// Set up progress property listeners
		this.progressProperty().addListener((observable, oldValue, newValue) -> {
			if (this.progressBar != null) {
				this.progressBar.setValue(newValue.intValue() * 100);
			}
		});

		// Set up handlers for task completion
		this.setOnSucceeded(event -> done());
		this.setOnFailed(event -> done());
		this.setOnCancelled(event -> {
			Logger LOGGER = getLogger();
			String taskName = getTaskName();
			LOGGER.info(taskName + ": cancelled");
			System.out.println(format.format(new Date()) + "\t" + taskName + ": task cancelled");
		});

		initial();
	}

	/**
	 * Initialize the task
	 */
	protected abstract void initial();

	/**
	 * Set the log builder
	 * 
	 * @param builder The configuration builder
	 */
	public void setLogBuilder(ConfigurationBuilder<BuiltConfiguration> builder) {
		this.builder = builder;
	}

	/**
	 * Get the task name
	 * 
	 * @return The task name
	 */
	protected abstract String getTaskName();

	/**
	 * Get the logger
	 * 
	 * @return The logger
	 */
	protected abstract Logger getLogger();

	/**
	 * This method is called when the task completes, whether successfully or not.
	 * It replaces the SwingWorker's done() method.
	 */
	protected void done() {
		Logger LOGGER = getLogger();
		String taskName = getTaskName();

		boolean finished = false;
		try {
			finished = get();
		} catch (InterruptedException | ExecutionException e) {
			LOGGER.error(taskName + ": failed", e);
			System.err.println(format.format(new Date()) + "\t" + taskName + ": task failed");
		}

		if (finished) {
			LOGGER.info(taskName + ": finished");
			System.out.println(format.format(new Date()) + "\t" + taskName + ": task finished");

			if (nextWork != null) {
				Thread thread = new Thread(nextWork);
				thread.setDaemon(true);
				thread.start();
			} else {
				updateProgress(1.0, 1.0);
			}
		} else {
			LOGGER.info(taskName + ": failed");
			System.out.println(format.format(new Date()) + "\t" + taskName + ": task failed");
		}
	}

	/**
	 * Compatibility method to start the task execution. Equivalent to SwingWorker's
	 * execute() method.
	 */
	public void execute() {
		Thread thread = new Thread(this);
		thread.setDaemon(true);
		thread.start();
	}

	/**
	 * Compatibility method to check if the task is done
	 * 
	 * @return true if the task is done, false otherwise
	 */
	public boolean isDone() {
		return this.getState() == State.SUCCEEDED || this.getState() == State.FAILED
				|| this.getState() == State.CANCELLED;
	}

	/**
	 * Set progress using the old 0-100 scale for compatibility with SwingWorker
	 * 
	 * @param progress value from 0 to 100
	 */
	protected void setProgress(int progress) {
		updateProgress(progress / 100.0, 1.0);
	}

	/**
	 * Update a message on the JavaFX UI thread
	 * 
	 * @param message The message to display
	 */
	protected void updateUIMessage(final String message) {
		Platform.runLater(() -> {
			updateMessage(message);
		});
	}
}
