package uk.ac.quadram.metapilot.gui.common;

import javafx.geometry.Insets;
import javafx.scene.control.*;
import javafx.scene.layout.*;

public class IOPanel extends VBox {
    public TableView<String> fileTable;
    public TextField microbiomeField;
    public TextField hostField;
    public CheckBox appendHostCheckBox;

    public IOPanel() {
        setPadding(new Insets(10));
        setSpacing(10);

        // Raw files section
        TitledPane rawFilesPane = new TitledPane("Raw Files", createRawFilesPanel());
        rawFilesPane.setExpanded(true);

        // Database section
        TitledPane databasePane = new TitledPane("Database", createDatabasePanel());
        databasePane.setExpanded(true);

        getChildren().addAll(rawFilesPane, databasePane);
    }

    private VBox createRawFilesPanel() {
        VBox panel = new VBox(5);
        panel.setPadding(new Insets(5));

        // File selection
        HBox fileSelectionBox = new HBox(10);
        Button addFilesButton = new Button("Add Files");
        Button addDirectoryButton = new Button("Add Directory");
        fileSelectionBox.getChildren().addAll(addFilesButton, addDirectoryButton);

        // File table
        fileTable = new TableView<>();
        TableColumn<String, String> fileNameColumn = new TableColumn<>("File Name");
        TableColumn<String, String> experimentNameColumn = new TableColumn<>("Experiment Name");
        fileTable.getColumns().addAll(fileNameColumn, experimentNameColumn);

        panel.getChildren().addAll(fileSelectionBox, fileTable);
        return panel;
    }

    private VBox createDatabasePanel() {
        VBox panel = new VBox(5);
        panel.setPadding(new Insets(5));

        // Microbiome database
        HBox microbiomeBox = new HBox(10);
        Label microbiomeLabel = new Label("Microbiome Database:");
        microbiomeField = new TextField();
        Button microbiomeBrowseButton = new Button("Browse");
        microbiomeBox.getChildren().addAll(microbiomeLabel, microbiomeField, microbiomeBrowseButton);

        // Host database
        HBox hostBox = new HBox(10);
        Label hostLabel = new Label("Host Database:");
        hostField = new TextField();
        Button hostBrowseButton = new Button("Browse");
        appendHostCheckBox = new CheckBox("Append Host Database");
        hostBox.getChildren().addAll(hostLabel, hostField, hostBrowseButton, appendHostCheckBox);

        panel.getChildren().addAll(microbiomeBox, hostBox);
        return panel;
    }
} 