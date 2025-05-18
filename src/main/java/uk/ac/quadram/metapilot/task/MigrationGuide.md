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

## Migration Steps

1. Change class declaration from extending `SwingWorker<Boolean, String>` to `Task<Boolean>`
2. Replace `doInBackground()` with `call()`
3. Replace `publish()` calls with `updateMessage()`
4. Replace `setProgress(0-100)` with `updateProgress(0.0-1.0, 1.0)`
5. Replace `cancel(true)` with `cancel()`
6. Update imports to use JavaFX classes

## Example

### Original SwingWorker code:

```java
public class MyTask extends SwingWorker<Boolean, String> {
    
    @Override
    protected Boolean doInBackground() throws Exception {
        setProgress(0);
        publish("Starting task...");
        
        // Do work
        
        setProgress(100);
        return true;
    }
    
    @Override
    protected void process(List<String> chunks) {
        // Update UI with the messages
    }
    
    @Override
    protected void done() {
        // Clean up after completion
    }
}
```

### Migrated to JavaFX Task:

```java
public class MyTask extends Task<Boolean> {
    
    @Override
    protected Boolean call() throws Exception {
        updateProgress(0, 100);
        updateMessage("Starting task...");
        
        // Do work
        
        updateProgress(100, 100);
        return true;
    }
    
    // For compatibility with old code
    public Boolean get() throws InterruptedException, ExecutionException {
        try {
            return this.getValue();
        } catch (Exception e) {
            throw new ExecutionException(e);
        }
    }
    
    // For compatibility with old code
    public void execute() {
        Thread thread = new Thread(this);
        thread.setDaemon(true);
        thread.start();
    }
}
```

## Handling UI Updates

In JavaFX, UI updates from background tasks are automatically handled through properties:

```java
Task<Boolean> task = new MyTask();
// Listen for state changes
task.setOnSucceeded(event -> {
    // Handle successful completion
});
task.setOnFailed(event -> {
    // Handle failure
});

// Bind progress to progress indicator/bar
progressBar.progressProperty().bind(task.progressProperty());
// Bind message to label
statusLabel.textProperty().bind(task.messageProperty());

// Start the task
new Thread(task).start();
// Or using the compatibility method:
// task.execute();
```

## Important Notes

1. JavaFX applications have a strong separation of UI and background threads. All UI updates must happen on the JavaFX Application Thread.
2. `Task` has more built-in support for properties that can be bound directly to UI controls.
3. When updating progress in JavaFX, you typically use values from 0.0 to 1.0, but we've maintained using 0-100 for compatibility.
4. The JavaFX Task's `updateProgress()` method expects a pair of numbers representing "progress" and "total". 