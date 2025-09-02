package uk.ac.quadram.metapilot.task.common.gui;

import javafx.collections.FXCollections;
import javafx.geometry.Insets;
import javafx.scene.control.*;
import javafx.scene.layout.*;
import javafx.scene.text.Font;
import javafx.scene.text.FontWeight;
import javafx.scene.paint.Color;

/**
 * Independent taxonomic analysis panel component that can be reused across different workflow panels.
 * Provides UI controls for configuring taxonomic analysis parameters including:
 * - Taxa selection (Eukaryota, Bacteria, Archaea, Viruses, Viroids, Unclassified, Other)
 * - Ignore rank selection
 * - Environmental samples filtering
 */
public class TaxonomicAnalysisPanel extends VBox {
    
    // Taxonomic analysis components
    private CheckBox eukaryotaCheckBox;
    private CheckBox bacteriaCheckBox;
    private CheckBox archaeaCheckBox;
    private CheckBox virusesCheckBox;
    private CheckBox viroidsCheckBox;
    private CheckBox unclassifiedCheckBox;
    private CheckBox otherCheckBox;
    private ComboBox<String> ignoreRankComboBox;
    private CheckBox environmentalSamplesCheckBox;
    
    // Callback interface for parameter updates
    public interface ParameterUpdateCallback {
        void onParameterUpdate();
    }
    
    private ParameterUpdateCallback updateCallback;
    
    /**
     * Constructor for TaxonomicAnalysisPanel
     */
    public TaxonomicAnalysisPanel() {
        this(null);
    }
    
    /**
     * Constructor for TaxonomicAnalysisPanel with callback
     * @param updateCallback Callback to be called when parameters are updated
     */
    public TaxonomicAnalysisPanel(ParameterUpdateCallback updateCallback) {
        this.updateCallback = updateCallback;
        initializeUI();
    }
    
    /**
     * Initialize the UI components
     */
    private void initializeUI() {
        setSpacing(8);
        setPadding(new Insets(8));
        setStyle("-fx-background-color: #f0fff0; -fx-background-radius: 5;");
        
        Label titleLabel = createSectionTitle("Taxonomic Analysis");
        
        // Taxa selection (horizontal compact grid - more columns, fewer rows)
        Label taxaLabel = new Label("Taxa from which top nodes showed in result:");
        taxaLabel.setStyle("-fx-font-weight: bold; -fx-font-size: 11px;");
        
        GridPane taxaGrid = new GridPane();
        taxaGrid.setHgap(8);
        taxaGrid.setVgap(5);
        
        // Initialize checkboxes
        eukaryotaCheckBox = new CheckBox("Eukaryota");
        bacteriaCheckBox = new CheckBox("Bacteria");
        archaeaCheckBox = new CheckBox("Archaea");
        virusesCheckBox = new CheckBox("Viruses");
        viroidsCheckBox = new CheckBox("Viroids");
        unclassifiedCheckBox = new CheckBox("Unclassified");
        otherCheckBox = new CheckBox("Other");
        
        // Add action listeners for parameter updates
        addParameterUpdateListener(eukaryotaCheckBox);
        addParameterUpdateListener(bacteriaCheckBox);
        addParameterUpdateListener(archaeaCheckBox);
        addParameterUpdateListener(virusesCheckBox);
        addParameterUpdateListener(viroidsCheckBox);
        addParameterUpdateListener(unclassifiedCheckBox);
        addParameterUpdateListener(otherCheckBox);
        
        taxaGrid.add(eukaryotaCheckBox, 0, 0);
        taxaGrid.add(bacteriaCheckBox, 1, 0);
        taxaGrid.add(archaeaCheckBox, 2, 0);
        taxaGrid.add(virusesCheckBox, 3, 0);
        taxaGrid.add(viroidsCheckBox, 4, 0);
        taxaGrid.add(unclassifiedCheckBox, 5, 0);
        taxaGrid.add(otherCheckBox, 6, 0);
        
        // Ignore rank - compact horizontal layout
        HBox rankBox = new HBox(8);
        rankBox.setAlignment(javafx.geometry.Pos.CENTER_LEFT);
        Label rankLabel = new Label("Ignore blanks below rank:");
        rankLabel.setStyle("-fx-font-weight: bold; -fx-font-size: 11px;");
        ignoreRankComboBox = new ComboBox<>(FXCollections.observableArrayList(
            "Superkingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"));
        ignoreRankComboBox.setValue("Species");
        ignoreRankComboBox.setMinWidth(80);
        ignoreRankComboBox.setPrefWidth(100);
        ignoreRankComboBox.setMaxWidth(Double.MAX_VALUE);
        ignoreRankComboBox.setOnAction(e -> notifyParameterUpdate());
        rankBox.getChildren().addAll(rankLabel, ignoreRankComboBox);
        
        environmentalSamplesCheckBox = new CheckBox("Ignore environmental samples");
        environmentalSamplesCheckBox.setOnAction(e -> notifyParameterUpdate());
        
        // Restore defaults button
        Button restoreDefaultsButton = new Button("Restore Defaults");
        restoreDefaultsButton.setOnAction(e -> restoreDefaults());
        
        getChildren().addAll(titleLabel, taxaLabel, taxaGrid, rankBox, environmentalSamplesCheckBox, restoreDefaultsButton);
    }
    
    /**
     * Add parameter update listener to a checkbox
     */
    private void addParameterUpdateListener(CheckBox checkBox) {
        checkBox.setOnAction(e -> notifyParameterUpdate());
    }
    
    /**
     * Notify callback of parameter update
     */
    private void notifyParameterUpdate() {
        if (updateCallback != null) {
            updateCallback.onParameterUpdate();
        }
    }
    
    /**
     * Create a section title label
     */
    private Label createSectionTitle(String title) {
        Label label = new Label(title);
        label.setFont(Font.font("System", FontWeight.BOLD, 14));
        label.setTextFill(Color.DARKBLUE);
        return label;
    }
    
    /**
     * Restore default taxonomic settings
     */
    public void restoreDefaults() {
        eukaryotaCheckBox.setSelected(true);
        bacteriaCheckBox.setSelected(true);
        archaeaCheckBox.setSelected(true);
        virusesCheckBox.setSelected(true);
        viroidsCheckBox.setSelected(false);
        unclassifiedCheckBox.setSelected(false);
        otherCheckBox.setSelected(false);
        ignoreRankComboBox.setValue("Family");
        environmentalSamplesCheckBox.setSelected(true);
        notifyParameterUpdate();
    }
    
    // Getters for taxonomic settings
    public boolean isEukaryota() { return eukaryotaCheckBox.isSelected(); }
    public boolean isBacteria() { return bacteriaCheckBox.isSelected(); }
    public boolean isArchaea() { return archaeaCheckBox.isSelected(); }
    public boolean isViruses() { return virusesCheckBox.isSelected(); }
    public boolean isViroids() { return viroidsCheckBox.isSelected(); }
    public boolean isUnclassified() { return unclassifiedCheckBox.isSelected(); }
    public boolean isOther() { return otherCheckBox.isSelected(); }
    public String getIgnoreRank() { return ignoreRankComboBox.getValue(); }
    public boolean isIgnoreEnvironmentalSamples() { return environmentalSamplesCheckBox.isSelected(); }
    
    // Setters for taxonomic settings
    public void setEukaryota(boolean value) { eukaryotaCheckBox.setSelected(value); }
    public void setBacteria(boolean value) { bacteriaCheckBox.setSelected(value); }
    public void setArchaea(boolean value) { archaeaCheckBox.setSelected(value); }
    public void setViruses(boolean value) { virusesCheckBox.setSelected(value); }
    public void setViroids(boolean value) { viroidsCheckBox.setSelected(value); }
    public void setUnclassified(boolean value) { unclassifiedCheckBox.setSelected(value); }
    public void setOther(boolean value) { otherCheckBox.setSelected(value); }
    public void setIgnoreRank(String value) { ignoreRankComboBox.setValue(value); }
    public void setIgnoreEnvironmentalSamples(boolean value) { environmentalSamplesCheckBox.setSelected(value); }
    
    /**
     * Set the parameter update callback
     */
    public void setParameterUpdateCallback(ParameterUpdateCallback callback) {
        this.updateCallback = callback;
    }
    
    /**
     * Get a summary of current taxonomic settings for display
     */
    public String getTaxonomicSummary() {
        StringBuilder sb = new StringBuilder();
        sb.append("=== TAXONOMIC ANALYSIS ===\n");
        sb.append("Eukaryota: ").append(isEukaryota() ? "Yes" : "No").append("\n");
        sb.append("Bacteria: ").append(isBacteria() ? "Yes" : "No").append("\n");
        sb.append("Archaea: ").append(isArchaea() ? "Yes" : "No").append("\n");
        sb.append("Viruses: ").append(isViruses() ? "Yes" : "No").append("\n");
        sb.append("Viroids: ").append(isViroids() ? "Yes" : "No").append("\n");
        sb.append("Unclassified: ").append(isUnclassified() ? "Yes" : "No").append("\n");
        sb.append("Other: ").append(isOther() ? "Yes" : "No").append("\n");
        sb.append("Ignore Rank: ").append(getIgnoreRank()).append("\n");
        sb.append("Ignore Environmental Samples: ").append(isIgnoreEnvironmentalSamples() ? "Yes" : "No").append("\n");
        return sb.toString();
    }
}
