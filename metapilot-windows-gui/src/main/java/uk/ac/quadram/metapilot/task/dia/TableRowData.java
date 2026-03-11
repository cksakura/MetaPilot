package uk.ac.quadram.metapilot.task.dia;

import javafx.beans.property.SimpleStringProperty;
import javafx.beans.property.StringProperty;

/**
 * Data class for JavaFX TableView rows in real-time DIA workflow
 */
public class TableRowData {
    private final StringProperty fileName;
    private final StringProperty status;
    private final StringProperty details;

    public TableRowData() {
        this.fileName = new SimpleStringProperty();
        this.status = new SimpleStringProperty();
        this.details = new SimpleStringProperty();
    }

    public TableRowData(String fileName, String status, String details) {
        this.fileName = new SimpleStringProperty(fileName);
        this.status = new SimpleStringProperty(status);
        this.details = new SimpleStringProperty(details);
    }

    // Getters for properties (for JavaFX binding)
    public StringProperty fileNameProperty() { return fileName; }
    public StringProperty statusProperty() { return status; }
    public StringProperty detailsProperty() { return details; }

    // Getters for values
    public String getFileName() { return fileName.get(); }
    public String getStatus() { return status.get(); }
    public String getDetails() { return details.get(); }

    // Setters
    public void setFileName(String fileName) { this.fileName.set(fileName); }
    public void setStatus(String status) { this.status.set(status); }
    public void setDetails(String details) { this.details.set(details); }

    // Helper method to find row by fileName
    public boolean matchesFileName(String fileName) {
        return this.fileName.get() != null && this.fileName.get().equals(fileName);
    }
}
