# MetaPilot Swing to JavaFX Migration Plan

## Current Status

The MetaPilot application is currently in a transition state with a hybrid architecture:

1. The main application entry point (`MetaPilotApp.java`) is already a JavaFX application.
2. The main application provides a modern JavaFX UI for the launcher/workflow selector.
3. When a specific workflow is launched, it creates a Swing-based JFrame (e.g., `MetaPilotMainFrameDia`).
4. Some components have already been migrated to JavaFX (e.g., `MagDbConfigPanelJavaFX`).

## Migration Approach

We will follow a modular, incremental approach to migrate the remaining Swing components to JavaFX:

1. **Create JavaFX equivalents** of the main frame classes first
2. **Convert individual panels** one by one
3. **Update workflow launching logic** in MetaPilotApp
4. **Test extensively** after each major component is migrated

## Core Components to Migrate

### 1. Main Frames

| Swing Class | JavaFX Equivalent |
|-------------|-------------------|
| MetaPilotMainFrameDia | MetaPilotMainFrameDiaFX |
| MetaPilotMainFrameMag | MetaPilotMainFrameMagFX |
| MetaPilotPFindMainFrame | MetaPilotPFindMainFrameFX |
| MetaPilotMainFrameDiaTg | MetaPilotMainFrameDiaTgFX |
| MetaPilotMainFrameDiaRT | MetaPilotMainFrameDiaRTFX |

### 2. Main Panels

| Swing Class | JavaFX Equivalent |
|-------------|-------------------|
| MetaPilotMainPanelDia | MetaPilotMainPanelDiaFX |
| MetaPilotMainPanelMag | MetaPilotMainPanelMagFX |
| MetaPilotMainPanelDiaRT | MetaPilotMainPanelDiaRTFX |

### 3. Dialog Classes

| Swing Class | JavaFX Equivalent |
|-------------|-------------------|
| MetaLabSourceDialogDia | MetaLabSourceDialogDiaFX |
| MetaLabSourcePathDialog | MetaLabSourcePathDialogFX |
| MetaPilotLicenseFrame | MetaPilotLicenseDialogFX |

### 4. Component Panels

| Swing Class | JavaFX Equivalent |
|-------------|-------------------|
| MagDbStatPanel | MagDbStatPanelFX |
| MagDbPepLibPanel | MagDbPepLibPanelFX |
| MagDbModelPanel | MagDbModelPanelFX |
| LibSearchPanel | LibSearchPanelFX |
| DiaNNSearchLibPanel | DiaNNSearchLibPanelFX |
| DiaNNCreateLibPanel | DiaNNCreateLibPanelFX |
| MetaLabIOPanel | MetaLabIOPanelFX |
| MetaLabMetaPanel | MetaLabMetaPanelFX |
| StatusPanel | StatusPanelFX |

## Technical Conversion Details

### JFrame → Stage

Convert JFrame instances to JavaFX Stage:

```java
// Swing
JFrame frame = new JFrame("Title");
frame.setSize(800, 600);
frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
frame.setVisible(true);

// JavaFX
Stage stage = new Stage();
stage.setTitle("Title");
stage.setWidth(800);
stage.setHeight(600);
stage.setOnCloseRequest(e -> Platform.exit());
stage.show();
```

### JPanel → Pane (or subclasses)

Convert JPanel to appropriate JavaFX containers:

```java
// Swing
JPanel panel = new JPanel();
panel.setLayout(new BorderLayout());

// JavaFX
BorderPane pane = new BorderPane();
```

### MigLayout → JavaFX Layouts

Replace MigLayout with appropriate JavaFX layouts:

```java
// Swing with MigLayout
panel.setLayout(new MigLayout("", "[grow]", "[grow]"));
panel.add(component, "cell 0 0");

// JavaFX 
GridPane gridPane = new GridPane();
gridPane.add(component, 0, 0);
// Or use other layouts like VBox, HBox, BorderPane as appropriate
```

### JButton → Button

```java
// Swing
JButton button = new JButton("Click Me");
button.addActionListener(e -> doSomething());

// JavaFX
Button button = new Button("Click Me");
button.setOnAction(e -> doSomething());
```

### SwingWorker → Task

Continue using the JavaFX Task for background processing (as you've already started):

```java
// Swing
SwingWorker<Boolean, Void> worker = new SwingWorker<>() {
    @Override
    protected Boolean doInBackground() throws Exception {
        // do work
        return true;
    }
    
    @Override
    protected void done() {
        // handle completion
    }
};
worker.execute();

// JavaFX
Task<Boolean> task = new Task<>() {
    @Override
    protected Boolean call() throws Exception {
        // do work
        return true;
    }
};
task.setOnSucceeded(e -> {
    // handle completion
});
new Thread(task).start();
```

### JProgressBar → ProgressBar

```java
// Swing
JProgressBar progressBar = new JProgressBar(0, 100);
progressBar.setValue(50);

// JavaFX
ProgressBar progressBar = new ProgressBar(0.5); // Takes values 0.0 to 1.0
```

### JTable → TableView

```java
// Swing
JTable table = new JTable(model);

// JavaFX
TableView<YourDataType> tableView = new TableView<>();
tableView.setItems(observableList);
tableView.getColumns().add(new TableColumn<>("Column Name"));
```

## Implementation Order

1. Start with the common utility components (status panels, progress indicators)
2. Convert the main panel classes next
3. Implement the main frame classes that will host these panels
4. Update the workflow launch logic in MetaPilotApp
5. Convert dialog classes
6. Finally address any remaining specialized components

## Testing Strategy

1. Create unit tests for each converted component
2. Implement integration tests for component interactions
3. Manually test the UI for appearance and behavior
4. Perform end-to-end workflow testing

## Resources

- [JavaFX API Documentation](https://openjfx.io/javadoc/21/)
- [JavaFX CSS Reference Guide](https://openjfx.io/javadoc/21/javafx.graphics/javafx/scene/doc-files/cssref.html)
- [Scene Builder](https://gluonhq.com/products/scene-builder/) for visual layout design 