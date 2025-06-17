package uk.ac.quadram.metapilot.gui.common;

import javafx.geometry.Insets;
import javafx.scene.control.*;
import javafx.scene.layout.*;

public class MetadataPanel extends VBox {
    public TableView<String> metadataTable;
    public Button updateButton;

    public MetadataPanel() {
        setPadding(new Insets(10));
        setSpacing(10);

        Label titleLabel = new Label("Metadata");
        titleLabel.setStyle("-fx-font-weight: bold;");

        metadataTable = new TableView<>();
        TableColumn<String, String> experimentColumn = new TableColumn<>("Experiment Name");
        TableColumn<String, String> sampleColumn = new TableColumn<>("Sample Info");
        metadataTable.getColumns().addAll(experimentColumn, sampleColumn);

        updateButton = new Button("Update");

        getChildren().addAll(titleLabel, metadataTable, updateButton);
    }
} 