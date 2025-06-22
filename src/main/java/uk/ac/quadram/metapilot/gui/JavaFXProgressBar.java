package uk.ac.quadram.metapilot.gui;

import javafx.application.Platform;
import javafx.geometry.Pos;
import javafx.scene.control.Label;
import javafx.scene.control.ProgressBar;
import javafx.scene.layout.StackPane;

/**
 * A wrapper for JavaFX ProgressBar that provides similar functionality to
 * Swing's JProgressBar to ease the transition between the two frameworks.
 * 
 * This class wraps a JavaFX ProgressBar and a Label to provide the string
 * painting functionality of JProgressBar.
 * 
 * @author Kai Cheng
 */
public class JavaFXProgressBar extends StackPane {
    
    private ProgressBar progressBar;
    private Label label;
    private boolean stringPainted = false;
    private int min = 0;
    private int max = 100;
    
    /**
     * Creates a new progress bar with initial value 0
     */
    public JavaFXProgressBar() {
        this(0, 100, 0);
    }
    
    /**
     * Creates a new progress bar with the specified min, max and initial value
     * 
     * @param min The minimum value (usually 0)
     * @param max The maximum value (usually 100)
     * @param initialValue The initial value
     */
    public JavaFXProgressBar(int min, int max, int initialValue) {
        this.min = min;
        this.max = max;
        
        // Create the progress bar
        progressBar = new ProgressBar(convertToProgress(initialValue));
        progressBar.setMaxWidth(Double.MAX_VALUE);
        
        // Create the label for string display
        label = new Label("");
        label.setAlignment(Pos.CENTER);
        label.setStyle("-fx-text-fill: black;");
        
        // Add both to the stack pane
        getChildren().addAll(progressBar, label);
    }
    
    /**
     * Convert the int value (0-100) to a progress value (0.0-1.0)
     * 
     * @param value The value to convert
     * @return The progress value
     */
    private double convertToProgress(int value) {
        return (double)(value - min) / (max - min);
    }
    
    /**
     * Set whether the progress bar should paint a string
     * 
     * @param b True if the progress bar should paint a string
     */
    public void setStringPainted(boolean b) {
        this.stringPainted = b;
        if (!stringPainted) {
            label.setText("");
        }
    }
    
    /**
     * Check if the progress bar is painting a string
     * 
     * @return True if the progress bar is painting a string
     */
    public boolean isStringPainted() {
        return stringPainted;
    }
    
    /**
     * Set the string to display on the progress bar
     * 
     * @param string The string to display
     */
    public void setString(String string) {
        if (stringPainted) {
            Platform.runLater(() -> label.setText(string));
        }
    }
    
    /**
     * Get the displayed string
     * 
     * @return The displayed string
     */
    public String getString() {
        return label.getText();
    }
    
    /**
     * Set the current value of the progress bar
     * 
     * @param value The value to set (0-100)
     */
    public void setValue(int value) {
        Platform.runLater(() -> {
            progressBar.setProgress(convertToProgress(value));
            
            // Update the percentage display if string is painted but no custom string is set
            if (stringPainted && label.getText().isEmpty()) {
                label.setText(value + "%");
            }
        });
    }
    
    /**
     * Get the current value of the progress bar
     * 
     * @return The current value (0-100)
     */
    public int getValue() {
        return (int) (progressBar.getProgress() * (max - min) + min);
    }
    
    /**
     * Set the minimum value of the progress bar
     * 
     * @param min The minimum value
     */
    public void setMinimum(int min) {
        this.min = min;
        setValue(getValue()); // Refresh display
    }
    
    /**
     * Get the minimum value of the progress bar
     * 
     * @return The minimum value
     */
    public int getMinimum() {
        return min;
    }
    
    /**
     * Set the maximum value of the progress bar
     * 
     * @param max The maximum value
     */
    public void setMaximum(int max) {
        this.max = max;
        setValue(getValue()); // Refresh display
    }
    
    /**
     * Get the maximum value of the progress bar
     * 
     * @return The maximum value
     */
    public int getMaximum() {
        return max;
    }
    
    /**
     * Set the progress bar to indeterminate mode
     * 
     * @param indeterminate True for indeterminate mode, false for determinate
     */
    public void setIndeterminate(boolean indeterminate) {
        Platform.runLater(() -> progressBar.setProgress(indeterminate ? ProgressBar.INDETERMINATE_PROGRESS : convertToProgress(getValue())));
    }
    
    /**
     * Check if the progress bar is in indeterminate mode
     * 
     * @return True if the progress bar is in indeterminate mode
     */
    public boolean isIndeterminate() {
        return progressBar.getProgress() == ProgressBar.INDETERMINATE_PROGRESS;
    }
    
    /**
     * Get the underlying JavaFX ProgressBar
     * 
     * @return The JavaFX ProgressBar
     */
    public ProgressBar getProgressBar() {
        return progressBar;
    }
} 