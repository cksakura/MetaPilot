package uk.ac.quadram.metapilot.gui.common;

import javafx.scene.control.TextArea;
import javafx.scene.layout.VBox;
import javafx.geometry.Insets;

public class ConsolePanel extends VBox {
    private TextArea consoleOutput;

    public ConsolePanel() {
        setPadding(new Insets(10));
        setSpacing(10);

        consoleOutput = new TextArea();
        consoleOutput.setEditable(false);
        consoleOutput.setWrapText(true);
        consoleOutput.setPrefRowCount(20);
        consoleOutput.setStyle("-fx-font-family: monospace;");

        getChildren().add(consoleOutput);
    }

    public void appendText(String text) {
        consoleOutput.appendText(text + "\n");
    }

    public void clear() {
        consoleOutput.clear();
    }
} 