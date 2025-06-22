package uk.ac.quadram.metapilot.task.common.gui;

import javafx.geometry.Insets;
import javafx.geometry.Pos;
import javafx.scene.control.*;
import javafx.scene.layout.*;
import javafx.animation.KeyFrame;
import javafx.animation.Timeline;
import javafx.util.Duration;
import java.lang.management.ManagementFactory;
import java.lang.management.MemoryMXBean;
import java.lang.management.MemoryUsage;

/**
 * Reusable console panel that can be used across different MetaPilot workflows.
 * Includes console output, progress bars, parameter view, and control buttons.
 */
public class ConsolePanel extends VBox {
    // UI Components
    private TextArea consoleOutput;
    private TextArea parameterView;
    private ProgressBar overallProgressBar;
    private ProgressBar currentTaskProgressBar;
    private Label currentTaskLabel;
    private Button startButton;
    private Button stopButton;
    
    // Status panel components
    private Label statusTimeLabel;
    private Label statusCpuLabel;
    private Label statusMemLabel;
    private Timeline statusTimeline;
    private long taskStartTimeMillis;
    
    // Callbacks for start/stop actions
    private Runnable startCallback;
    private Runnable stopCallback;
    
    public ConsolePanel() {
        initializeUI();
    }
    
    private void initializeUI() {
        setSpacing(10);
        setPadding(new Insets(10));
        
        // Parameter View
        VBox paramViewBox = new VBox(5);
        paramViewBox.getChildren().add(new Label("Parameter View:"));
        parameterView = new TextArea();
        parameterView.setEditable(false);
        parameterView.setPrefHeight(200);
        parameterView.setStyle("-fx-font-family: monospace;");
        VBox.setVgrow(parameterView, Priority.ALWAYS);
        paramViewBox.getChildren().add(parameterView);
        VBox.setVgrow(paramViewBox, Priority.ALWAYS);
        
        // Progress Bars
        VBox progressBox = new VBox(5);
        currentTaskLabel = new Label("Current Task: None");
        overallProgressBar = new ProgressBar(0);
        overallProgressBar.setMaxWidth(Double.MAX_VALUE);
        currentTaskProgressBar = new ProgressBar(0);
        currentTaskProgressBar.setMaxWidth(Double.MAX_VALUE);
        progressBox.getChildren().addAll(
            currentTaskLabel,
            new Label("Overall Progress:"),
            overallProgressBar,
            new Label("Current Task Progress:"),
            currentTaskProgressBar
        );
        
        // Console Output
        VBox consoleBox = new VBox(5);
        consoleBox.getChildren().add(new Label("Console Output:"));
        consoleOutput = new TextArea();
        consoleOutput.setEditable(false);
        consoleOutput.setPrefHeight(250);
        consoleOutput.setStyle("-fx-font-family: monospace;");
        VBox.setVgrow(consoleOutput, Priority.ALWAYS);
        consoleBox.getChildren().add(consoleOutput);
        VBox.setVgrow(consoleBox, Priority.ALWAYS);
        
        // Control Buttons
        HBox controlBox = new HBox(10);
        startButton = new Button("Start Analysis");
        stopButton = new Button("Stop");
        stopButton.setDisable(true);
        startButton.setOnAction(e -> handleStart());
        stopButton.setOnAction(e -> handleStop());
        controlBox.getChildren().addAll(startButton, stopButton);
        
        // Status Panel (bottom)
        HBox statusPanel = new HBox(20);
        statusPanel.setPadding(new Insets(5, 0, 0, 0));
        statusPanel.setAlignment(Pos.CENTER_LEFT);
        statusTimeLabel = new Label("Time: 00:00:00");
        statusCpuLabel = new Label("CPU: 0.0%");
        statusMemLabel = new Label("Memory: 0 / 0 MB");
        statusPanel.getChildren().addAll(statusTimeLabel, statusCpuLabel, statusMemLabel);
        
        // Layout
        getChildren().addAll(
            paramViewBox,
            progressBox,
            consoleBox,
            controlBox,
            statusPanel
        );
        VBox.setVgrow(paramViewBox, Priority.ALWAYS);
        VBox.setVgrow(consoleBox, Priority.ALWAYS);
        VBox.setVgrow(this, Priority.ALWAYS);
    }
    
    // Public methods for external control
    public void setStartCallback(Runnable callback) {
        this.startCallback = callback;
    }
    
    public void setStopCallback(Runnable callback) {
        this.stopCallback = callback;
    }
    
    public void setParameterText(String text) {
        if (parameterView != null) {
            parameterView.setText(text);
        }
    }
    
    public void appendConsoleText(String text) {
        if (consoleOutput != null) {
            consoleOutput.appendText(text);
        }
    }
    
    public void clearConsole() {
        if (consoleOutput != null) {
            consoleOutput.clear();
        }
    }
    
    public void setCurrentTaskLabel(String text) {
        if (currentTaskLabel != null) {
            currentTaskLabel.setText(text);
        }
    }
    
    public ProgressBar getOverallProgressBar() {
        return overallProgressBar;
    }
    
    public ProgressBar getCurrentTaskProgressBar() {
        return currentTaskProgressBar;
    }
    
    public void startStatusTimer() {
        taskStartTimeMillis = System.currentTimeMillis();
        if (statusTimeline != null) statusTimeline.stop();
        statusTimeline = new Timeline(new KeyFrame(Duration.seconds(1), e -> updateStatusPanel()));
        statusTimeline.setCycleCount(Timeline.INDEFINITE);
        statusTimeline.play();
    }
    
    public void stopStatusTimer() {
        if (statusTimeline != null) {
            statusTimeline.stop();
        }
    }
    
    public void setButtonsEnabled(boolean startEnabled, boolean stopEnabled) {
        if (startButton != null) startButton.setDisable(!startEnabled);
        if (stopButton != null) stopButton.setDisable(!stopEnabled);
    }
    
    private void handleStart() {
        if (startCallback != null) {
            setButtonsEnabled(false, true);
            clearConsole();
            startStatusTimer();
            startCallback.run();
        }
    }
    
    private void handleStop() {
        if (stopCallback != null) {
            setButtonsEnabled(true, false);
            stopStatusTimer();
            stopCallback.run();
        }
    }
    
    private void updateStatusPanel() {
        if (statusTimeLabel != null && taskStartTimeMillis > 0) {
            long elapsed = System.currentTimeMillis() - taskStartTimeMillis;
            long hours = elapsed / 3600000;
            long minutes = (elapsed % 3600000) / 60000;
            long seconds = (elapsed % 60000) / 1000;
            statusTimeLabel.setText(String.format("Time: %02d:%02d:%02d", hours, minutes, seconds));
        }
        
        if (statusCpuLabel != null) {
            // Simplified CPU usage display
            statusCpuLabel.setText("CPU: N/A");
        }
        
        if (statusMemLabel != null) {
            MemoryMXBean memoryBean = ManagementFactory.getMemoryMXBean();
            MemoryUsage heapUsage = memoryBean.getHeapMemoryUsage();
            long usedMemory = heapUsage.getUsed() / (1024 * 1024); // Convert to MB
            statusMemLabel.setText("Memory: " + usedMemory + " MB");
        }
    }
} 