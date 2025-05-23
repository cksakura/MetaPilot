# Migration Guide: SwingWorker to JavaFX Task

This document guides you through the process of migrating task classes from Swing's `SwingWorker` to JavaFX's `Task`.

## Key Differences

| SwingWorker | JavaFX Task | Notes |
|-------------|-------------|-------|
| `doInBackground()` | `call()` | Main method where background work is performed |
| `process()` | N/A | No direct equivalent - use `updateMessage()` instead |
| `publish()` | `updateMessage()` | Send progress messages |
| `setProgress(int)` | `updateProgress(double, double)` | Update progress (JavaFX uses values from 0.0 to 1.0) |
| `execute()` | N/A | Create Thread and call start() (added as a compatibility method) |
| `get()` | `getValue()` | Get the result (added as a compatibility method) |
| `cancel(boolean)` | `cancel()` | Cancel the task |
| `done()` | `setOnSucceeded()`, `setOnFailed()`, `setOnCancelled()` | Called after task completion |

## Progress Bar Migration

We've created a `JavaFXProgressBar` class that wraps the JavaFX `ProgressBar` control to provide the same API as Swing's `JProgressBar`:

```java
// Old Swing code
JProgressBar progressBar = new JProgressBar(0, 100);
progressBar.setStringPainted(true);
progressBar.setValue(50);
progressBar.setString("Processing...");

// New JavaFX code
JavaFXProgressBar progressBar = new JavaFXProgressBar(0, 100, 0);
progressBar.setStringPainted(true);
progressBar.setValue(50);
progressBar.setString("Processing...");
```

## Migration Steps

1. Change class declaration from extending `SwingWorker<Boolean, String>` to `MetaAbstractTask`
2. Replace `doInBackground()` with `call()`
3. Replace `publish()` calls with `updateMessage()` or `updateUIMessage()`
4. Either keep using `setProgress(0-100)` (compatibility method) or use `updateProgress(0.0-1.0, 1.0)`
5. Replace `cancel(true)` with `cancel()`
6. Update imports to use JavaFX classes
7. Replace `JProgressBar` with `JavaFXProgressBar`
8. Keep using the `done()` method (preserved for compatibility) or use event handlers

## Base Class Overview

The MetaAbstractTask extends JavaFX Task and provides compatibility functionality:

- **Progress Tracking**: Both the JavaFX progress model (0.0-1.0) and SwingWorker model (0-100)
- **Task Lifecycle**: Overrides for `call()`, `done()`, and automatic event handling
- **Progress Reporting**: Updates progress bars and status messages
- **Compatibility Methods**: `get()`, `execute()`, `isDone()`, and `setProgress(int)`
- **UI Thread Safety**: `updateUIMessage()` for thread-safe UI updates

## Example Migration

### Original SwingWorker code:

```java
public class MyTask extends SwingWorker<Boolean, String> {
    
    private JProgressBar progressBar;
    
    public MyTask(JProgressBar progressBar) {
        this.progressBar = progressBar;
        addPropertyChangeListener(evt -> {
            if ("progress".equals(evt.getPropertyName())) {
                progressBar.setValue((Integer) evt.getNewValue());
            }
        });
    }
    
    @Override
    protected Boolean doInBackground() throws Exception {
        setProgress(0);
        publish("Starting task...");
        
        // Do some work
        for (int i = 0; i < 100; i++) {
            Thread.sleep(100);
            setProgress(i + 1);
        }
        
        return true;
    }
    
    @Override
    protected void process(List<String> chunks) {
        for (String chunk : chunks) {
            System.out.println(chunk);
        }
    }
    
    @Override
    protected void done() {
        try {
            boolean result = get();
            if (result) {
                System.out.println("Task completed successfully");
            } else {
                System.out.println("Task failed");
            }
        } catch (Exception e) {
            e.printStackTrace();
        }
    }
}
```

### Migrated to JavaFX Task:

```java
public class MyTask extends MetaAbstractTask {
    
    private JavaFXProgressBar progressBar;
    
    public MyTask(JavaFXProgressBar progressBar) {
        super(null, null, progressBar, null, null);
        // Progress bar binding is handled by the superclass constructor
    }
    
    @Override
    protected void initial() {
        // Any initialization code
    }
    
    @Override
    protected String getTaskName() {
        return "My Task";
    }
    
    @Override
    protected Logger getLogger() {
        return LogManager.getLogger(MyTask.class);
    }
    
    @Override
    protected Boolean call() throws Exception {
        setProgress(0);  // Can keep using setProgress(int) for compatibility
        updateMessage("Starting task...");
        
        // Do some work
        for (int i = 0; i < 100; i++) {
            Thread.sleep(100);
            setProgress(i + 1);  // Works the same as in SwingWorker
        }
        
        return true;
    }
    
    // done() method is inherited from MetaAbstractTask
    // If you want custom behavior, you can override it
}
```

## Using with FXML

When working with FXML-defined UIs, you can bind the progress property directly to progress bars:

```java
@FXML
private ProgressBar progressBar;

private void setupTask() {
    MyTask task = new MyTask();
    progressBar.progressProperty().bind(task.progressProperty());
    
    // For message updates
    task.messageProperty().addListener((obs, oldValue, newValue) -> {
        statusLabel.setText(newValue);
    });
    
    new Thread(task).start();
}
```

## Common Issues to Watch For

1. **Task Execution**: Use the `execute()` compatibility method or manually create a thread and start it
2. **Progress Updates**: You can continue using `setProgress(int)` (0-100 range) or use JavaFX's `updateProgress(double, double)` (0.0-1.0 range)
3. **UI Thread Updates**: Use `updateUIMessage()` for Platform.runLater updates or bind to message property
4. **Task Cancellation**: Use `cancel()` method (no mayInterruptIfRunning parameter)
5. **Task Results**: Use `get()` compatibility method to retrieve results
6. **Exception Handling**: Override `done()` or use event handlers

## Tips for Testing

1. Create simple unit tests for migrated tasks
2. Test task execution, cancellation, and progress reporting
3. Verify that UI updates correctly with task progress
4. Ensure proper lifecycle handling (task completion, failure, cancellation)
5. Test task chaining to ensure subsequent tasks execute properly 