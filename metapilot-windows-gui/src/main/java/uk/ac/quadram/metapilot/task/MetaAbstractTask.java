package uk.ac.quadram.metapilot.task;

import java.text.SimpleDateFormat;
import java.util.Date;
import java.util.concurrent.ExecutionException;

import javafx.application.Platform;
import javafx.concurrent.Task;
import javafx.concurrent.Worker;
import javafx.scene.control.ProgressBar;

import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.core.config.builder.api.ConfigurationBuilder;
import org.apache.logging.log4j.core.config.builder.impl.BuiltConfiguration;

import uk.ac.quadram.metapilot.resources.MetaResourceParameters;
import uk.ac.quadram.metapilot.task.par.MetaBaseParameter;

/**
 * Abstract Base Class for MetaPilot Background Tasks
 * 
 * <h2>Purpose</h2> This abstract class serves as the foundation for all
 * computational tasks in MetaPilot, providing a unified framework for
 * background processing, progress tracking, and lifecycle management while
 * maintaining compatibility with the original SwingWorker design.
 * 
 * <h2>Architecture</h2>
 * <ul>
 * <li><strong>JavaFX Task Integration:</strong> Extends JavaFX
 * Task&lt;Boolean&gt; for modern asynchronous processing</li>
 * <li><strong>Progress Tracking:</strong> Built-in progress monitoring with UI
 * updates</li>
 * <li><strong>Lifecycle Management:</strong> Standardized task initialization,
 * execution, and completion handling</li>
 * <li><strong>Error Handling:</strong> Comprehensive exception management with
 * logging</li>
 * <li><strong>Chain Execution:</strong> Support for sequential task
 * execution</li>
 * </ul>
 * 
 * <h2>Task Categories</h2>
 * <p>
 * This base class supports various task types in MetaPilot:
 * </p>
 * <ul>
 * <li><strong>Database Search Tasks:</strong> Mass spectrometry database
 * searches</li>
 * <li><strong>Analysis Tasks:</strong> Protein identification and
 * quantification</li>
 * <li><strong>Annotation Tasks:</strong> Taxonomic and functional
 * annotation</li>
 * <li><strong>Data Processing Tasks:</strong> File conversion and
 * preprocessing</li>
 * <li><strong>Database Management Tasks:</strong> Reference database
 * operations</li>
 * </ul>
 * 
 * <h2>Implementation Pattern</h2>
 * <p>
 * Subclasses must implement:
 * </p>
 * <ul>
 * <li>{@link #initial()} - Task-specific initialization</li>
 * <li>{@link #call()} - Main task execution logic</li>
 * <li>{@link #getTaskName()} - Descriptive task name for logging</li>
 * <li>{@link #getLogger()} - Task-specific logger instance</li>
 * </ul>
 * 
 * <h2>Progress Reporting</h2>
 * <p>
 * Tasks can report progress through:
 * </p>
 * <ul>
 * <li>{@link #updateProgress(double, double)} - JavaFX progress updates</li>
 * <li>{@link #setProgress(int)} - Legacy 0-100 scale progress</li>
 * <li>{@link #updateMessage(String)} - Status message updates</li>
 * </ul>
 * 
 * <h2>Threading Model</h2>
 * <ul>
 * <li><strong>Background Execution:</strong> All tasks run in daemon
 * threads</li>
 * <li><strong>UI Thread Safety:</strong> Progress updates are marshaled to
 * JavaFX Application Thread</li>
 * <li><strong>Resource Management:</strong> Automatic cleanup on task
 * completion</li>
 * </ul>
 * 
 * <h2>Task Chaining</h2>
 * <p>
 * Tasks can be chained for sequential execution:
 * </p>
 * 
 * <pre>{@code
 *   MetaAbstractTask searchTask = new DatabaseSearchTask(...);
 *   MetaAbstractTask annotationTask = new AnnotationTask(...);
 *   searchTask.setNextWork(annotationTask);
 *   searchTask.execute();
 * }</pre>
 * 
 * <h2>Configuration</h2>
 * <p>
 * Tasks are configured with:
 * </p>
 * <ul>
 * <li><strong>metaPar:</strong> Analysis parameters specific to the
 * workflow</li>
 * <li><strong>mrp:</strong> Resource paths for external tools and
 * databases</li>
 * <li><strong>progressBar:</strong> UI component for progress
 * visualization</li>
 * </ul>
 * 
 * @author Kai Cheng - Quadram Institute Bioscience
 * @version 2.0
 * @since 1.0
 * 
 * @see javafx.concurrent.Task
 * @see uk.ac.quadram.metapilot.task.par.MetaBaseParameter
 * @see uk.ac.quadram.metapilot.resources.MetaResourceParameters
 * @see uk.ac.quadram.metapilot.gui.JavaFXProgressBar
 */
public abstract class MetaAbstractTask extends Task<Boolean> {

	protected MetaBaseParameter metaPar;
	protected MetaResourceParameters mrp;
	protected ProgressBar progressBar;
	protected Logger logger;
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
	public MetaAbstractTask(MetaBaseParameter metaPar, MetaResourceParameters mrp, ProgressBar progressBar,
			Logger logger) {
		this.metaPar = metaPar;
		this.mrp = mrp;
		this.progressBar = progressBar;
		this.logger = logger;
		// Set up progress property listeners
		this.progressProperty().addListener((observable, oldValue, newValue) -> {
			if (this.progressBar != null) {
				this.progressBar.setProgress(newValue.doubleValue());
			}
		});

		// Set up handlers for task completion
		this.setOnSucceeded(event -> done());
		this.setOnFailed(event -> done());
		this.setOnCancelled(event -> {
			String taskName = getTaskName();
			logger.info("[{}] task cancelled", taskName);
		});
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
	public MetaAbstractTask(MetaBaseParameter metaPar, MetaResourceParameters mrp, ProgressBar progressBar,
			Logger logger, Task<Boolean> nextWork) {
		this.metaPar = metaPar;
		this.mrp = mrp;
		this.progressBar = progressBar;
		this.logger = logger;
		this.nextWork = nextWork;

		// Set up progress property listeners
		this.progressProperty().addListener((observable, oldValue, newValue) -> {
			if (this.progressBar != null) {
				this.progressBar.setProgress(newValue.doubleValue());
			}
		});

		// Set up handlers for task completion
		this.setOnSucceeded(event -> done());
		this.setOnFailed(event -> done());
		this.setOnCancelled(event -> {
			String taskName = getTaskName();
			logger.info("[{}] task cancelled", taskName);
		});
	}

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
	protected Logger getLogger() {
		return logger;
	}

	/**
	 * This method is called when the task completes, whether successfully or not.
	 * It replaces the SwingWorker's done() method.
	 */
	protected void done() {
		String taskName = getTaskName();
		boolean finished = false;
		try {
			finished = get();
		} catch (InterruptedException | ExecutionException e) {
			logger.error(taskName + ": failed", e);
			System.err.println(format.format(new Date()) + "\t" + taskName + ": task failed");
		}

		if (finished) {
			logger.info("[{}] task finished", taskName);

			if (nextWork != null) {
				Thread thread = new Thread(nextWork);
				thread.setDaemon(true);
				thread.start();
			} else {
				updateProgress(1.0, 1.0);
			}
		} else {
			logger.info("[{}] task failed", taskName);
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
		return this.getState() == Worker.State.SUCCEEDED || this.getState() == Worker.State.FAILED
				|| this.getState() == Worker.State.CANCELLED;
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
