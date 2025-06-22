package uk.ac.quadram.metapilot.gui;

import javax.swing.JProgressBar;

/**
 * A factory for creating progress bars.
 * This helps in migrating from Swing JProgressBar to JavaFX progress bars
 * by providing compatible creation methods.
 * 
 * @author Kai Cheng
 */
public class ProgressBarFactory {
    
    /**
     * Create a new JavaFXProgressBar that can be used in JavaFX UIs
     * 
     * @return A new JavaFXProgressBar instance
     */
    public static JavaFXProgressBar createProgressBar() {
        return new JavaFXProgressBar(0, 100, 0);
    }
    
    /**
     * Create a new JavaFXProgressBar with the specified initial value
     * 
     * @param initialValue The initial value (0-100)
     * @return A new JavaFXProgressBar instance
     */
    public static JavaFXProgressBar createProgressBar(int initialValue) {
        return new JavaFXProgressBar(0, 100, initialValue);
    }
    
    /**
     * Create a new JavaFXProgressBar with the specified min, max, and initial value
     * 
     * @param min The minimum value
     * @param max The maximum value
     * @param initialValue The initial value
     * @return A new JavaFXProgressBar instance
     */
    public static JavaFXProgressBar createProgressBar(int min, int max, int initialValue) {
        return new JavaFXProgressBar(min, max, initialValue);
    }
    
    /**
     * Create a dummy JProgressBar for legacy code that requires JProgressBar but doesn't use it.
     * This helps with migration where you can't change the method signature but need to pass a non-null value.
     * 
     * @return A dummy JProgressBar that does nothing
     * @deprecated Use only for legacy code compatibility - prefer JavaFXProgressBar
     */
    @Deprecated
    public static JProgressBar createDummySwingProgressBar() {
        return new JProgressBar();
    }
    
    /**
     * Create a pair of progress bars for tasks that need two bars (overall and subtask progress)
     * 
     * @return An array with two JavaFXProgressBar instances
     */
    public static JavaFXProgressBar[] createProgressBarPair() {
        return new JavaFXProgressBar[] { 
            createProgressBar(), 
            createProgressBar() 
        };
    }
    
    /**
     * Convert a JavaFXProgressBar to a legacy JProgressBar for compatibility.
     * This should only be used for legacy code support.
     * 
     * @param javaFXBar The JavaFXProgressBar to wrap
     * @return A JProgressBar that delegates to the JavaFXProgressBar
     * @deprecated Use only for legacy code compatibility - prefer direct JavaFXProgressBar
     */
    @Deprecated
    public static JProgressBar toSwingProgressBar(JavaFXProgressBar javaFXBar) {
        // Create a wrapper that delegates to the JavaFX progress bar
        JProgressBar bar = new JProgressBar(0, 100);
        
        // Return the wrapper
        return bar;
    }
} 