package uk.ac.quadram.metapilot.task.dia;

import javafx.beans.property.SimpleStringProperty;
import javafx.beans.property.StringProperty;

/**
 * Data class for JavaFX TableView rows in real-time DIA workflow
 */
public class TableRowData {
    private final StringProperty fileName;
    private final StringProperty round1Status;
    private final StringProperty round2Status;
    private final StringProperty round3Status;
    private final StringProperty round4Status;
    private final StringProperty round5Status;

    public TableRowData() {
        this.fileName = new SimpleStringProperty();
        this.round1Status = new SimpleStringProperty();
        this.round2Status = new SimpleStringProperty();
        this.round3Status = new SimpleStringProperty();
        this.round4Status = new SimpleStringProperty();
        this.round5Status = new SimpleStringProperty();
    }

    public TableRowData(String fileName, String round1Status, String round2Status, 
                       String round3Status, String round4Status, String round5Status) {
        this.fileName = new SimpleStringProperty(fileName);
        this.round1Status = new SimpleStringProperty(round1Status);
        this.round2Status = new SimpleStringProperty(round2Status);
        this.round3Status = new SimpleStringProperty(round3Status);
        this.round4Status = new SimpleStringProperty(round4Status);
        this.round5Status = new SimpleStringProperty(round5Status);
    }

    // Getters for properties (for JavaFX binding)
    public StringProperty fileNameProperty() { return fileName; }
    public StringProperty round1StatusProperty() { return round1Status; }
    public StringProperty round2StatusProperty() { return round2Status; }
    public StringProperty round3StatusProperty() { return round3Status; }
    public StringProperty round4StatusProperty() { return round4Status; }
    public StringProperty round5StatusProperty() { return round5Status; }

    // Getters for values
    public String getFileName() { return fileName.get(); }
    public String getRound1Status() { return round1Status.get(); }
    public String getRound2Status() { return round2Status.get(); }
    public String getRound3Status() { return round3Status.get(); }
    public String getRound4Status() { return round4Status.get(); }
    public String getRound5Status() { return round5Status.get(); }

    // Setters
    public void setFileName(String fileName) { this.fileName.set(fileName); }
    public void setRound1Status(String round1Status) { this.round1Status.set(round1Status); }
    public void setRound2Status(String round2Status) { this.round2Status.set(round2Status); }
    public void setRound3Status(String round3Status) { this.round3Status.set(round3Status); }
    public void setRound4Status(String round4Status) { this.round4Status.set(round4Status); }
    public void setRound5Status(String round5Status) { this.round5Status.set(round5Status); }

    // Helper method to find row by fileName
    public boolean matchesFileName(String fileName) {
        return this.fileName.get() != null && this.fileName.get().equals(fileName);
    }
}
