package uk.ac.quadram.metapilot.task.common.gui;

import javafx.geometry.Insets;
import javafx.scene.control.*;
import javafx.scene.layout.*;
import javafx.stage.FileChooser;
import javafx.collections.FXCollections;
import javafx.collections.ObservableList;
import javafx.scene.text.Font;
import javafx.scene.text.FontWeight;
import javafx.scene.paint.Color;
import uk.ac.quadram.metapilot.task.mag.MagDbItem;
import uk.ac.quadram.metapilot.task.mag.MagDbItem.MagVersion;
import uk.ac.quadram.metapilot.task.mag.MagDbConfigIO;

import java.io.File;
import java.util.ArrayList;
import java.util.List;

/**
 * Reusable MAG Database Settings Panel
 * 
 * This panel provides a consistent interface for selecting MAG databases,
 * versions, and host database options across different workflow panels.
 */
public class MagDatabaseSettingsPanel extends VBox {
    
    // UI Components
    private ComboBox<MagDbItem> magDatabaseComboBox;
    private ComboBox<MagVersion> magVersionComboBox;
    private TextField hostDBField;
    private CheckBox appendHostDBCheckBox;
    private Button browseHostDBButton;
    
    // Data
    private ObservableList<MagDbItem> availableMagDatabases = FXCollections.observableArrayList();
    
    // Getter for available databases
    public ObservableList<MagDbItem> getAvailableMagDatabases() {
        return availableMagDatabases;
    }
    
    // Callback for parameter updates
    private Runnable parameterUpdateCallback;
    
    // Whether to show host database options (DDA workflows) or hide them (DIA workflows)
    private boolean showHostDatabaseOptions;
    
    public MagDatabaseSettingsPanel() {
        this(true); // Default to showing host database options for backward compatibility
    }
    
    public MagDatabaseSettingsPanel(boolean showHostDatabaseOptions) {
        this.showHostDatabaseOptions = showHostDatabaseOptions;
        initializeUI();
    }
    
    private void initializeUI() {
        setSpacing(20);
        setPadding(new Insets(10));
        setStyle("-fx-background-color: #f8f8f8; -fx-border-color: #cccccc; -fx-border-radius: 5; -fx-border-width: 1;");
        
        // Title
        Label titleLabel = createSectionTitle("MAG Database Settings");
        
        // MAG Database ComboBox
        HBox magDBBox = new HBox(10);
        magDatabaseComboBox = new ComboBox<>(availableMagDatabases);
        magDatabaseComboBox.setPromptText("Select MAG database");
        magDatabaseComboBox.setPrefWidth(150);
        magDatabaseComboBox.setOnAction(e -> {
            updateAvailableVersions();
            if (parameterUpdateCallback != null) {
                parameterUpdateCallback.run();
            }
        });
        magDBBox.getChildren().addAll(new Label("MAG Database:"), magDatabaseComboBox);
        
        // MAG Version ComboBox
        magVersionComboBox = new ComboBox<>();
        magVersionComboBox.setPromptText("Select MAG version");
        magVersionComboBox.setPrefWidth(100);
        magVersionComboBox.setDisable(true); // Initially disabled, will be enabled after database selection
        magVersionComboBox.setOnAction(e -> {
            if (parameterUpdateCallback != null) {
                parameterUpdateCallback.run();
            }
        });
        magDBBox.getChildren().addAll(new Label("MAG Version:"), magVersionComboBox);
        
        // Host Database Options (only show for DDA workflows)
        if (showHostDatabaseOptions) {
            appendHostDBCheckBox = new CheckBox("Append Host Database");
            appendHostDBCheckBox.setOnAction(e -> {
                updateHostDBVisibility();
                if (parameterUpdateCallback != null) {
                    parameterUpdateCallback.run();
                }
            });
            magDBBox.getChildren().add(appendHostDBCheckBox);
            
            hostDBField = new TextField();
            hostDBField.setPromptText("Select host database file");
            hostDBField.setPrefWidth(600);
            hostDBField.setDisable(true);
            browseHostDBButton = new Button("Browse");
            browseHostDBButton.setDisable(true);
            browseHostDBButton.setOnAction(e -> browseHostDatabase());
            magDBBox.getChildren().addAll(new Label("Host Database:"), hostDBField, browseHostDBButton);
        } else {
            // For DIA workflows, initialize as null to avoid null pointer exceptions
            appendHostDBCheckBox = null;
            hostDBField = null;
            browseHostDBButton = null;
        }
        
        getChildren().addAll(titleLabel, magDBBox);
    }
    
    private Label createSectionTitle(String title) {
        Label label = new Label(title);
        label.setFont(Font.font("System", FontWeight.BOLD, 14));
        label.setTextFill(Color.DARKBLUE);
        return label;
    }
    
    private void updateHostDBVisibility() {
        if (appendHostDBCheckBox != null) {
            boolean selected = appendHostDBCheckBox.isSelected();
            if (hostDBField != null) {
                hostDBField.setDisable(!selected);
            }
            if (browseHostDBButton != null) {
                browseHostDBButton.setDisable(!selected);
            }
        }
    }
    
    private void browseHostDatabase() {
        if (hostDBField == null) return; // Safety check for DIA workflows
        
        FileChooser chooser = new FileChooser();
        chooser.setTitle("Select Host Database File");
        chooser.getExtensionFilters().add(
            new FileChooser.ExtensionFilter("FASTA Files", "*.fasta", "*.fa")
        );
        File file = chooser.showOpenDialog(getScene().getWindow());
        if (file != null) {
            hostDBField.setText(file.getAbsolutePath());
            if (parameterUpdateCallback != null) {
                parameterUpdateCallback.run();
            }
        }
    }
    
    /**
     * Set available MAG databases
     */
    public void setAvailableMagDatabases(MagDbItem[] magDbItems) {
        availableMagDatabases.clear();
        if (magDbItems != null) {
            for (MagDbItem item : magDbItems) {
                availableMagDatabases.add(item);
            }
        }
        
        // Set default selection if available
        if (!availableMagDatabases.isEmpty()) {
            magDatabaseComboBox.setValue(availableMagDatabases.get(0));
            updateAvailableVersions();
        }
    }
    
    /**
     * Update available versions when a MAG database is selected
     */
    private void updateAvailableVersions() {
        MagDbItem selectedDb = magDatabaseComboBox.getValue();
        
        magVersionComboBox.getItems().clear();
        magVersionComboBox.setDisable(true);
        
        if (selectedDb != null && selectedDb.getVersionList() != null) {
            for (MagVersion version : selectedDb.getVersionList()) {
                magVersionComboBox.getItems().add(version);
            }
            
            if (!magVersionComboBox.getItems().isEmpty()) {
                magVersionComboBox.setDisable(false);
                // Set default version (first available)
                magVersionComboBox.setValue(magVersionComboBox.getItems().get(0));
            }
        }
    }
    
    /**
     * Load MAG databases from configuration files
     */
    public MagDbItem[] loadMagDatabasesFromConfig() {
        try {
            // Try to load from the default MAG database configuration file
            File configFile = new File("MAGDB_1_0.json");
            if (configFile.exists()) {
                MagDbItem[][] magDbItems = MagDbConfigIO.parse(configFile);
                if (magDbItems != null && magDbItems.length > 0) {
                    // Combine MGnify and local databases
                    List<MagDbItem> allDatabases = new ArrayList<>();
                    for (MagDbItem[] dbArray : magDbItems) {
                        if (dbArray != null) {
                            for (MagDbItem db : dbArray) {
                                if (db != null) {
                                    // Only include databases that have local files available
                                    if (db.getVersionList() != null) {
                                        for (MagVersion version : db.getVersionList()) {
                                            if (version.isLocalAvailable()) {
                                                allDatabases.add(db);
                                                break;
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                    return allDatabases.toArray(new MagDbItem[0]);
                }
            }
        } catch (Exception e) {
            System.err.println("ERROR: Failed to load MAG databases from config: " + e.getMessage());
            e.printStackTrace();
        }
        return new MagDbItem[0];
    }
    
    // Getters for current selections
    public MagDbItem getSelectedMagDatabase() {
        MagVersion selectedVersion = magVersionComboBox.getValue();
        if (selectedVersion != null) {
            MagDbItem magDbItem = magDatabaseComboBox.getValue();
            if (magDbItem != null) {
                magDbItem.setUsedVersion(
                    selectedVersion
                );
                return magDbItem;
            }
        }
        return null;
    }
    
    public MagVersion getSelectedMagVersion() {
        return magVersionComboBox.getValue();
    }
    
    public boolean isAppendHostDb() {
        return appendHostDBCheckBox != null ? appendHostDBCheckBox.isSelected() : false;
    }
    
    public String getHostDatabase() {
        return hostDBField != null ? hostDBField.getText() : "";
    }
    
    // Setters for programmatic updates
    public void setSelectedMagDatabase(MagDbItem magDb) {
        magDatabaseComboBox.setValue(magDb);
        updateAvailableVersions();
    }
    
    public void setSelectedMagVersion(MagVersion version) {
        magVersionComboBox.setValue(version);
    }
    
    public void setAppendHostDb(boolean append) {
        if (appendHostDBCheckBox != null) {
            appendHostDBCheckBox.setSelected(append);
            updateHostDBVisibility();
        }
    }
    
    public void setHostDatabase(String hostDb) {
        if (hostDBField != null) {
            hostDBField.setText(hostDb);
        }
    }
    
    // Callback setter
    public void setParameterUpdateCallback(Runnable callback) {
        this.parameterUpdateCallback = callback;
    }
    
    /**
     * Validate that all required fields are selected
     */
    public boolean validateSelection() {
        if (getSelectedMagDatabase() == null) {
            return false;
        }
        if (getSelectedMagVersion() == null) {
            return false;
        }
        // Only validate host database if host database options are shown (DDA workflows)
        if (showHostDatabaseOptions && isAppendHostDb() && (getHostDatabase() == null || getHostDatabase().trim().isEmpty())) {
            return false;
        }
        return true;
    }
    
    /**
     * Get validation error message
     */
    public String getValidationErrorMessage() {
        if (getSelectedMagDatabase() == null) {
            return "Please select a MAG database.";
        }
        if (getSelectedMagVersion() == null) {
            return "Please select a MAG version.";
        }
        // Only validate host database if host database options are shown (DDA workflows)
        if (showHostDatabaseOptions && isAppendHostDb() && (getHostDatabase() == null || getHostDatabase().trim().isEmpty())) {
            return "Please select a host database file.";
        }
        return null;
    }
}
