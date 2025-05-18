package uk.ac.quadram.metapilot.task;

import java.io.File;
import java.text.SimpleDateFormat;
import java.util.Date;
import java.util.concurrent.ExecutionException;

import javax.swing.JProgressBar;
import javafx.concurrent.Task;
import javafx.concurrent.WorkerStateEvent;
import javafx.event.EventHandler;

import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.core.config.builder.api.ConfigurationBuilder;
import org.apache.logging.log4j.core.config.builder.impl.BuiltConfiguration;

import uk.ac.quadram.metapilot.task.par.MetaParameter;
import uk.ac.quadram.metapilot.task.par.MetaSources;

/**
 * An abstract class of the tasks in MetaLab, which provides the progress bars, and handles
 * task completion using JavaFX's Task framework.
 * 
 * This class replaces MetaAbstractTask as part of the migration from Swing to JavaFX.
 * 
 * @author MEDTECHU0188
 *
 */
public abstract class JavaFXMetaAbstractTask extends Task<Boolean> {

    protected MetaParameter metaPar;
    protected MetaSources msv;
    protected JProgressBar bar1, bar2;
    protected Task<Boolean> nextWork;
    protected SimpleDateFormat format = new SimpleDateFormat("MM_dd_yyyy_HH_mm_ss");
    protected ConfigurationBuilder<BuiltConfiguration> builder;

    public JavaFXMetaAbstractTask(MetaParameter metaPar, MetaSources msv, JProgressBar bar1, JProgressBar bar2,
            Task<Boolean> nextWork) {
        this.metaPar = metaPar;
        this.msv = msv;
        this.bar1 = bar1;
        this.bar2 = bar2;
        this.nextWork = nextWork;

        // Set up progress property listeners
        this.progressProperty().addListener((observable, oldValue, newValue) -> {
            if (this.bar1 != null) {
                this.bar1.setValue(newValue.intValue() * 100);
            }
        });
        
        // Set up handlers for task completion
        this.setOnSucceeded(event -> handleTaskCompletion(event));
        this.setOnFailed(event -> handleTaskCompletion(event));
        this.setOnCancelled(event -> handleTaskCancelled(event));
        
        initial();
    }
    
    /**
     * Initialize the task
     */
    protected abstract void initial();
    
    /**
     * Set the log builder
     * @param builder The configuration builder
     */
    public void setLogBuilder(ConfigurationBuilder<BuiltConfiguration> builder) {
        this.builder = builder;
    }

    /**
     * Get the task name
     * @return The task name
     */
    protected abstract String getTaskName();

    /**
     * Get the logger
     * @return The logger
     */
    protected abstract Logger getLogger();

    /**
     * Handle task completion, whether successful or failed
     * @param event The worker state event
     */
    protected void handleTaskCompletion(WorkerStateEvent event) {
        Logger LOGGER = getLogger();
        String taskName = getTaskName();

        boolean finished = false;
        
        try {
            if (event.getEventType() == WorkerStateEvent.WORKER_STATE_SUCCEEDED) {
                finished = getValue();
            }
        } catch (Exception e) {
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
                updateProgress(100, 100);
            }
        } else {
            LOGGER.info(taskName + ": failed");
            System.out.println(format.format(new Date()) + "\t" + taskName + ": task failed");
        }
    }
    
    /**
     * Handle task cancellation
     * @param event The worker state event
     */
    protected void handleTaskCancelled(WorkerStateEvent event) {
        Logger LOGGER = getLogger();
        String taskName = getTaskName();
        
        LOGGER.info(taskName + ": cancelled");
        System.out.println(format.format(new Date()) + "\t" + taskName + ": task cancelled");
    }
    
    /**
     * Compatibility method to replace SwingWorker's get()
     * @return The result of the task execution
     * @throws InterruptedException if the task is interrupted
     * @throws ExecutionException if the task throws an exception
     */
    public Boolean get() throws InterruptedException, ExecutionException {
        try {
            // Wait for the task to complete and return the result
            return this.getValue();
        } catch (Exception e) {
            throw new ExecutionException(e);
        }
    }
    
    /**
     * Compatibility method to start the task execution.
     * Equivalent to SwingWorker's execute() method.
     */
    public void execute() {
        Thread thread = new Thread(this);
        thread.setDaemon(true);
        thread.start();
    }
    
    /**
     * Compatibility method to check if the task is done
     * @return true if the task is done, false otherwise
     */
    public boolean isDone() {
        return this.getState() == State.SUCCEEDED || 
               this.getState() == State.FAILED || 
               this.getState() == State.CANCELLED;
    }
} 