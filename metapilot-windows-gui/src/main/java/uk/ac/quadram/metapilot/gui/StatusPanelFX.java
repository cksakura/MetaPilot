package uk.ac.quadram.metapilot.gui;

import java.lang.management.ManagementFactory;
import java.lang.management.MemoryMXBean;
import java.text.DecimalFormat;

import javax.management.AttributeList;
import javax.management.InstanceNotFoundException;
import javax.management.MBeanServer;
import javax.management.MalformedObjectNameException;
import javax.management.ObjectName;
import javax.management.ReflectionException;

import javafx.animation.Animation;
import javafx.animation.KeyFrame;
import javafx.animation.Timeline;
import javafx.application.Application;
import javafx.geometry.Insets;
import javafx.geometry.Pos;
import javafx.scene.Scene;
import javafx.scene.control.Label;
import javafx.scene.layout.BorderPane;
import javafx.scene.layout.GridPane;
import javafx.scene.layout.HBox;
import javafx.scene.layout.Priority;
import javafx.scene.layout.Region;
import javafx.scene.paint.Color;
import javafx.stage.Stage;
import javafx.util.Duration;
import uk.ac.quadram.metapilot.core.tools.FormatTool;

/**
 * JavaFX implementation of the StatusPanel
 * Displays runtime metrics including elapsed time, CPU usage, and memory usage
 * 
 * @author Kai Cheng
 *
 */
public class StatusPanelFX extends GridPane {

    private DecimalFormat dfTime = new DecimalFormat("00");
    private DecimalFormat dfper = FormatTool.getDFP2();

    private Timeline timeline;
    private int counter = 0; // the duration

    private Label lblTime;
    private Label labelCpuUsed;
    private Label labelMemoryUsed;

    /**
     * Create the panel.
     */
    public StatusPanelFX(String version) {
        // Setup basic layout and styling
        this.setPadding(new Insets(10));
        this.setHgap(10);
        this.setVgap(5);
        this.getStyleClass().add("status-panel");
        
        // Create all labels
        Label lblRunTime = new Label("Running time:");
        lblTime = new Label(" ");
        
        Label lblCpu = new Label("CPU state:");
        labelCpuUsed = new Label(" ");
        
        Label lblMemory = new Label("Memory (used/max):");
        labelMemoryUsed = new Label(" ");
        
        Label lblVersion = new Label(version);
        lblVersion.setTextFill(Color.GRAY);
        
        // Add a spacer region that will push the version label to the right
        Region spacer = new Region();
        HBox.setHgrow(spacer, Priority.ALWAYS);
        
        // Add components to the grid
        // Row 0
        this.add(lblRunTime, 0, 0);
        this.add(lblTime, 1, 0);
        this.add(lblCpu, 2, 0);
        this.add(labelCpuUsed, 3, 0);
        this.add(lblMemory, 4, 0);
        this.add(labelMemoryUsed, 5, 0);
        this.add(spacer, 6, 0);
        this.add(lblVersion, 7, 0);
        
        // Column constraints to match the Swing layout
        for (int i = 0; i < 6; i += 2) {
            this.getColumnConstraints().add(new javafx.scene.layout.ColumnConstraints(100, 100, 100));
            this.getColumnConstraints().add(new javafx.scene.layout.ColumnConstraints(150, 200, 250));
        }
        
        // Spacer column that grows
        javafx.scene.layout.ColumnConstraints growingCol = new javafx.scene.layout.ColumnConstraints(50, 150, Double.MAX_VALUE);
        growingCol.setHgrow(Priority.ALWAYS);
        this.getColumnConstraints().add(growingCol);
        
        // Version column
        this.getColumnConstraints().add(new javafx.scene.layout.ColumnConstraints(100, 100, 100));
    }

    /**
     * Start the status monitoring timeline
     */
    public void start() {
        // Create a timeline that updates every second
        timeline = new Timeline(new KeyFrame(Duration.seconds(1), event -> {
            // Update time display
            int h = counter / 3600;
            int m = (counter - h * 3600) / 60;
            int s = counter - h * 3600 - m * 60;

            lblTime.setText(dfTime.format(h) + ":" + dfTime.format(m) + ":" + dfTime.format(s));
            counter++;

            // Update memory usage
            MemoryMXBean memoryMXBean = ManagementFactory.getMemoryMXBean();
            labelMemoryUsed.setText((memoryMXBean.getHeapMemoryUsage().getUsed() / 1024 / 1024) + " MB/"
                    + memoryMXBean.getHeapMemoryUsage().getMax() / 1024 / 1024 + " MB");

            // Update CPU usage
            MBeanServer mbs = ManagementFactory.getPlatformMBeanServer();
            ObjectName name = null;
            AttributeList list = null;
            
            try {
                name = ObjectName.getInstance("java.lang:type=OperatingSystem");
                list = mbs.getAttributes(name, new String[] { "SystemCpuLoad" });
            } catch (MalformedObjectNameException | NullPointerException | InstanceNotFoundException
                    | ReflectionException e) {
                e.printStackTrace();
            }

            if (list != null && !list.isEmpty()) {
                javax.management.Attribute att = (javax.management.Attribute) list.get(0);
                Double value = (Double) att.getValue();
                if (value > 0) {
                    labelCpuUsed.setText(dfper.format(value));
                }
            }
        }));
        
        timeline.setCycleCount(Animation.INDEFINITE);
        timeline.play();
    }

    /**
     * Stop the status monitoring timeline
     */
    public void stop() {
        if (timeline != null) {
            timeline.stop();
            timeline = null;
        }
    }

    /**
     * Main method for testing the component
     */
    public static void main(String[] args) {
        Application.launch(StatusPanelTestApp.class, args);
    }

    /**
     * Test application class
     */
    public static class StatusPanelTestApp extends Application {
        @Override
        public void start(Stage stage) {
            BorderPane root = new BorderPane();
            StatusPanelFX statusPanel = new StatusPanelFX("2.0.0");
            root.setCenter(statusPanel);
            
            Scene scene = new Scene(root, 800, 100);
            stage.setTitle("StatusPanelFX Test");
            stage.setScene(scene);
            stage.show();
            
            statusPanel.start();
        }
    }
} 