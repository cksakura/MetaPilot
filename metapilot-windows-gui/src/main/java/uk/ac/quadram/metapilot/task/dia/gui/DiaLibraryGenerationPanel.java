package uk.ac.quadram.metapilot.task.dia.gui;

import javafx.application.Platform;
import javafx.geometry.Insets;
import javafx.geometry.Pos;
import javafx.scene.control.*;
import javafx.scene.layout.*;
import javafx.collections.FXCollections;
import javafx.collections.ObservableList;
import javafx.scene.text.Font;
import javafx.scene.text.FontWeight;
import javafx.scene.paint.Color;
import javafx.scene.control.cell.PropertyValueFactory;
import javafx.beans.property.SimpleStringProperty;
import javafx.beans.property.StringProperty;
import javafx.beans.property.SimpleBooleanProperty;
import javafx.beans.property.BooleanProperty;
import javafx.concurrent.Task;
import uk.ac.quadram.metapilot.task.mag.MagDbItem;
import uk.ac.quadram.metapilot.task.mag.MagDbItem.MagPepLib;
import java.io.File;
import java.util.Optional;

import org.apache.logging.log4j.LogManager;

import javafx.animation.KeyFrame;
import javafx.animation.Timeline;
import javafx.util.Duration;
import uk.ac.quadram.metapilot.task.dia.DiaLibCreateTask;
import uk.ac.quadram.metapilot.task.common.gui.AbstractWorkflowPanel;
import uk.ac.quadram.metapilot.task.dia.par.MetaParameterDia;
import uk.ac.quadram.metapilot.resources.MetaResourceParameters;
import uk.ac.quadram.metapilot.task.mag.MagDbConfigIO;
import java.util.ArrayList;

public class DiaLibraryGenerationPanel extends VBox {

	// Library table components
	private TableView<LibraryRow> libraryTable;
	private ObservableList<LibraryRow> libraryData = FXCollections.observableArrayList();

	// Generation components
	private ToggleGroup inputTypeGroup;
	private TextField inputFileField;
	private ComboBox<String> enzymeComboBox;
	private ComboBox<Integer> missCleavageComboBox;
	private TextField libraryNameField;
	private VBox enzymeBox;

	// Current MAG database reference
	private MagDbItem currentMagDbItem;

	// Console output for feedback
	private TextArea consoleOutput;

	// Internal progress tracking
	private ProgressBar libraryProgressBar;
	private Label progressLabel;

	// Step tracking for clear progress indication
	private int currentStep = 0;
	private int totalSteps = 4;

	/** When true, sections are not added to this panel; parent uses getControlsSection() and getTableSection(). */
	private boolean splitMode;

	// Callback interface for parameter updates
	public interface ParameterUpdateCallback {
		void onParametersChanged();
	}

	private ParameterUpdateCallback parameterCallback;

	// Required parameters for DiaLibCreateTask
	private MetaParameterDia metaParameter;
	private MetaResourceParameters resourceParameters;
	private File configFile;
	private AbstractWorkflowPanel<MetaParameterDia> parentPanel;

	public DiaLibraryGenerationPanel() {
		initializeUI();
	}

	public DiaLibraryGenerationPanel(TextArea consoleOutput) {
		this.consoleOutput = consoleOutput;
		initializeUI();
	}

	/** For use in split layout (e.g. real-time panel): controls and table are placed in parent's columns. */
	public DiaLibraryGenerationPanel(TextArea consoleOutput, boolean splitMode) {
		this.consoleOutput = consoleOutput;
		this.splitMode = splitMode;
		initializeUI();
	}

	public DiaLibraryGenerationPanel(MetaParameterDia metaParameter, MetaResourceParameters resourceParameters,
			File configFile) {
		this.metaParameter = metaParameter;
		this.resourceParameters = resourceParameters;
		this.configFile = configFile;

		// CRITICAL: Initialize MagDbItem from metaParameter immediately
		if (metaParameter != null && metaParameter.getMagDbItem() != null) {
			this.currentMagDbItem = metaParameter.getMagDbItem();
			System.out.println("DEBUG: DiaLibraryGenerationPanel constructor - initialized MagDbItem: "
					+ currentMagDbItem.getCatalogueID());
		} else {
			System.err.println("ERROR: DiaLibraryGenerationPanel constructor - MetaParameterDia or MagDbItem is null!");
		}

		initializeUI();

		// Load available libraries after UI is initialized
		loadAvailableLibraries();
	}

	public DiaLibraryGenerationPanel(MetaParameterDia metaParameter, MetaResourceParameters resourceParameters,
			File configFile, TextArea consoleOutput) {
		this.metaParameter = metaParameter;
		this.resourceParameters = resourceParameters;
		this.configFile = configFile;
		this.consoleOutput = consoleOutput;

		// CRITICAL: Initialize MagDbItem from metaParameter immediately
		if (metaParameter != null && metaParameter.getMagDbItem() != null) {
			this.currentMagDbItem = metaParameter.getMagDbItem();
			System.out.println("DEBUG: DiaLibraryGenerationPanel constructor - initialized MagDbItem: "
					+ currentMagDbItem.getCatalogueID());
		} else {
			System.err.println("ERROR: DiaLibraryGenerationPanel constructor - MetaParameterDia or MagDbItem is null!");
		}

		initializeUI();

		// Load available libraries after UI is initialized
		loadAvailableLibraries();
	}

	public void setParameterUpdateCallback(ParameterUpdateCallback callback) {
		this.parameterCallback = callback;
	}

	public void setConsoleOutput(TextArea consoleOutput) {
		this.consoleOutput = consoleOutput;
	}

	public void setMagDbItem(MagDbItem magDbItem) {
		this.currentMagDbItem = magDbItem;
		if (magDbItem != null) {
			magDbItem.registerMarker4SpectralLibraryIfPresent();
		}
		loadAvailableLibraries();
	}

	public void setMetaParameter(MetaParameterDia metaParameter) {
		this.metaParameter = metaParameter;
	}

	public void setResourceParameters(MetaResourceParameters resourceParameters) {
		this.resourceParameters = resourceParameters;
	}

	public void setConfigFile(File configFile) {
		this.configFile = configFile;
	}

	public void setParentPanel(AbstractWorkflowPanel<MetaParameterDia> parentPanel) {
		this.parentPanel = parentPanel;
	}

	public void updateMetaParameter(MetaParameterDia metaParameter) {
		this.metaParameter = metaParameter;
	}

	private VBox controlsSectionVBox;
	private VBox tableSectionVBox;

	private void initializeUI() {
		setSpacing(10);
		setPadding(new Insets(10));
		setStyle("-fx-background-color: #fff8dc; -fx-background-radius: 5;");
		setMinHeight(400);
		VBox.setVgrow(this, Priority.ALWAYS);

		Label titleLabel = createSectionTitle("Spectral Library Generation");

		// Input type selection
		HBox inputTypeBox = createInputTypeSelection();

		// File selection
		HBox fileBox = createFileSelection();

		// Enzyme parameters (for protein database option)
		enzymeBox = createEnzymeParameters();

		// Generation controls
		HBox generationBox = createGenerationControls();

		// Progress and console section
		VBox progressSection = createProgressSection();

		// Available libraries table
		tableSectionVBox = createLibraryTable();

		// Left section: title + controls + progress
		controlsSectionVBox = new VBox(10);
		controlsSectionVBox.setPadding(new Insets(0));
		controlsSectionVBox.getChildren().addAll(titleLabel, inputTypeBox, fileBox, enzymeBox, generationBox,
				progressSection);
		VBox.setVgrow(progressSection, Priority.ALWAYS);

		if (!splitMode) {
			getChildren().addAll(controlsSectionVBox, tableSectionVBox);
		}
	}

	public VBox getControlsSection() {
		return controlsSectionVBox;
	}

	public VBox getTableSection() {
		return tableSectionVBox;
	}

	private HBox createInputTypeSelection() {
		HBox inputTypeBox = new HBox(8);
		Label inputLabel = new Label("Input Type for Library Generation:");
		inputLabel.setStyle("-fx-font-weight: bold;");

		inputTypeGroup = new ToggleGroup();
		RadioButton libraryRadio = new RadioButton("Existing Library File");
		RadioButton peptideListRadio = new RadioButton("Peptide List File");
		RadioButton proteinDbRadio = new RadioButton("Protein Database File");

		libraryRadio.setToggleGroup(inputTypeGroup);
		peptideListRadio.setToggleGroup(inputTypeGroup);
		proteinDbRadio.setToggleGroup(inputTypeGroup);
		libraryRadio.setSelected(true); // Default selection

		// Enable/disable enzyme parameters based on input type
		inputTypeGroup.selectedToggleProperty().addListener((obs, oldToggle, newToggle) -> {
			boolean enableEnzyme = newToggle == proteinDbRadio;
			enzymeComboBox.setDisable(!enableEnzyme);
			missCleavageComboBox.setDisable(!enableEnzyme);
		});

		inputTypeBox.getChildren().addAll(inputLabel, libraryRadio, peptideListRadio, proteinDbRadio);
		return inputTypeBox;
	}

	private HBox createFileSelection() {
		HBox fileBox = new HBox(10);
		fileBox.setAlignment(Pos.CENTER_LEFT);
		Label fileLabel = new Label("Input File:");
		fileLabel.setStyle("-fx-font-weight: bold;");
		inputFileField = new TextField();
		inputFileField.setPromptText("Select input file");
		inputFileField.setPrefWidth(250);
		inputFileField.setEditable(false);
		Button browseFileButton = new Button("Browse");
		browseFileButton.setOnAction(e -> browseInputFile());
		fileBox.getChildren().addAll(fileLabel, inputFileField, browseFileButton);
		return fileBox;
	}

	private VBox createEnzymeParameters() {
		VBox enzymeBox = new VBox(8);
		enzymeBox.setStyle("-fx-background-color: #f0f0f0; -fx-padding: 10; -fx-background-radius: 3;");
		Label enzymeTitle = new Label("Enzyme Parameters (for Protein Database):");
		enzymeTitle.setStyle("-fx-font-weight: bold;");

		GridPane enzymeGrid = new GridPane();
		enzymeGrid.setHgap(10);
		enzymeGrid.setVgap(5);

		// Enzyme name
		Label enzymeNameLabel = new Label("Enzyme:");
		enzymeComboBox = new ComboBox<>(
				FXCollections.observableArrayList("Trypsin/P", "Trypsin", "Lys-C", "Chymotrypsin", "AspN", "GluC"));
		enzymeComboBox.setValue("Trypsin/P");
		enzymeComboBox.setPrefWidth(120);

		// Miss cleavage
		Label missCleavageLabel = new Label("Max Miss Cleavages:");
		missCleavageComboBox = new ComboBox<>();
		missCleavageComboBox.getItems().addAll(0, 1, 2);
		missCleavageComboBox.setValue(1);
		missCleavageComboBox.setPrefWidth(80);

		enzymeGrid.add(enzymeNameLabel, 0, 0);
		enzymeGrid.add(enzymeComboBox, 1, 0);
		enzymeGrid.add(missCleavageLabel, 0, 1);
		enzymeGrid.add(missCleavageComboBox, 1, 1);

		enzymeBox.getChildren().addAll(enzymeTitle, enzymeGrid);

		// Enzyme parameters are always visible but initially disabled
		enzymeBox.setVisible(true);
		enzymeBox.setManaged(true);

		// Initially disable enzyme parameters since "Existing Library File" is selected by default
		enzymeComboBox.setDisable(true);
		missCleavageComboBox.setDisable(true);

		return enzymeBox;
	}

	private HBox createGenerationControls() {
		HBox generationBox = new HBox(10);
		generationBox.setAlignment(Pos.CENTER_LEFT);
		libraryNameField = new TextField();
		libraryNameField.setPromptText("Library name");
		libraryNameField.setPrefWidth(150);

		Button generateLibraryButton = new Button("Generate Library");
		generateLibraryButton.setPrefWidth(130);
		generateLibraryButton.setOnAction(e -> generateLibrary());

		generationBox.getChildren().addAll(new Label("Library Name:"), libraryNameField, generateLibraryButton);
		return generationBox;
	}

	private VBox createLibraryTable() {
		VBox tableSection = new VBox(10);

		// Header with controls
		HBox headerBox = new HBox(10);
		headerBox.setAlignment(Pos.CENTER_LEFT);
		Label availableLabel = new Label("Available Libraries (select for identification):");
		availableLabel.setStyle("-fx-font-weight: bold;");

		Button checkAvailableButton = new Button("Check Available");
		checkAvailableButton.setTooltip(new Tooltip("Auto-detect libraries in MAG database folder"));
		checkAvailableButton.setOnAction(e -> checkAvailableLibraries());

		Button deleteButton = new Button("Delete Selected");
		deleteButton.setOnAction(e -> deleteSelectedLibrary());

		Button testSelectButton = new Button("Test Select First");
		testSelectButton.setOnAction(e -> testSelectFirstLibrary());

		headerBox.getChildren().addAll(availableLabel, checkAvailableButton, deleteButton, testSelectButton);

		// Table
		libraryTable = new TableView<>(libraryData);
		libraryTable.setMinHeight(120);
		libraryTable.setPrefHeight(150);
		libraryTable.setMaxHeight(Double.MAX_VALUE);
		libraryTable.setEditable(true);
		libraryTable.setColumnResizePolicy(TableView.UNCONSTRAINED_RESIZE_POLICY);
		VBox.setVgrow(libraryTable, Priority.ALWAYS);

		// Use column (checkbox)
		TableColumn<LibraryRow, Boolean> useCol = new TableColumn<>("Use");
		useCol.setCellValueFactory(new PropertyValueFactory<>("use"));
		useCol.setCellFactory(column -> {
			return new TableCell<LibraryRow, Boolean>() {
				private final CheckBox checkBox = new CheckBox();

				{
					checkBox.setOnAction(e -> {
						LibraryRow row = getTableRow().getItem();
						if (row != null) {
							boolean newValue = checkBox.isSelected();
							System.out.println("DEBUG: Checkbox clicked - " + row.getLibName() + " -> " + newValue);
							row.setUse(newValue);
							notifyParameterUpdate();
						}
					});
				}

				@Override
				protected void updateItem(Boolean item, boolean empty) {
					super.updateItem(item, empty);
					if (empty || getTableRow() == null || getTableRow().getItem() == null) {
						setGraphic(null);
					} else {
						LibraryRow row = getTableRow().getItem();
						checkBox.setSelected(row.getUse());
						setGraphic(checkBox);
					}
				}
			};
		});
		useCol.setMinWidth(40);
		useCol.setPrefWidth(50);

		// Library Name column
		TableColumn<LibraryRow, String> nameCol = new TableColumn<>("Library Name");
		nameCol.setCellValueFactory(new PropertyValueFactory<>("libName"));
		nameCol.setMinWidth(120);
		nameCol.setPrefWidth(150);

		// Peptide Count column
		TableColumn<LibraryRow, String> pepCountCol = new TableColumn<>("Peptide Count");
		pepCountCol.setCellValueFactory(new PropertyValueFactory<>("pepCount"));
		pepCountCol.setMinWidth(80);
		pepCountCol.setPrefWidth(100);

		// Protease column
		TableColumn<LibraryRow, String> proteaseCol = new TableColumn<>("Protease");
		proteaseCol.setCellValueFactory(new PropertyValueFactory<>("protease"));
		proteaseCol.setMinWidth(80);
		proteaseCol.setPrefWidth(100);

		// Miss Cleavages column
		TableColumn<LibraryRow, String> missCol = new TableColumn<>("Miss Cleavages");
		missCol.setCellValueFactory(new PropertyValueFactory<>("miss"));
		missCol.setMinWidth(80);
		missCol.setPrefWidth(100);

		// Library type column (Microbiome vs Standalone)
		TableColumn<LibraryRow, String> typeCol = new TableColumn<>("Type");
		typeCol.setCellValueFactory(new PropertyValueFactory<>("libType"));
		typeCol.setMinWidth(100);
		typeCol.setPrefWidth(140);
		typeCol.setCellFactory(column -> {
			return new TableCell<LibraryRow, String>() {
				private final ComboBox<String> combo = new ComboBox<>(
						FXCollections.observableArrayList("Microbiome", "Standalone"));

				@Override
				protected void updateItem(String item, boolean empty) {
					super.updateItem(item, empty);
					if (empty || getTableRow() == null || getTableRow().getItem() == null) {
						setGraphic(null);
					} else {
						LibraryRow row = getTableRow().getItem();
						combo.setValue(row.getLibType());
						combo.setOnAction(e -> {
							row.setLibType(combo.getValue());
							if (row.getMagPepLib() != null) {
								row.getMagPepLib().setMatchToMag("Microbiome".equals(combo.getValue()));
							}
							notifyParameterUpdate();
						});
						setGraphic(combo);
					}
				}
			};
		});

		@SuppressWarnings("unchecked")
		TableColumn<LibraryRow, ?>[] columns = new TableColumn[] { useCol, nameCol, pepCountCol, proteaseCol, missCol, typeCol };
		libraryTable.getColumns().addAll(columns);

		tableSection.getChildren().addAll(headerBox, libraryTable);
		return tableSection;
	}

	private VBox createProgressSection() {
		VBox progressSection = new VBox(8);
		progressSection.setPadding(new Insets(10));
		progressSection.setStyle("-fx-background-color: #f8f8f8; -fx-background-radius: 5;");
		progressSection.setMinHeight(200); // Increased minimum height
		progressSection.setPrefHeight(400); // Set preferred height to make it bigger
		VBox.setVgrow(progressSection, Priority.ALWAYS); // Allow it to grow and fill space

		libraryProgressBar = new ProgressBar(0);
		libraryProgressBar.setMaxWidth(Double.MAX_VALUE);
		libraryProgressBar.setPrefHeight(25); // Made progress bar taller
		libraryProgressBar.setStyle("-fx-accent: #4CAF50;"); // Green color for progress

		// Progress label to show current step
		progressLabel = new Label("Ready to generate library");
		progressLabel.setStyle("-fx-font-weight: bold; -fx-font-size: 12px; -fx-text-fill: #333333;");
		progressLabel.setMaxWidth(Double.MAX_VALUE);

		Label consoleTitle = new Label("Console Output:");
		consoleTitle.setStyle("-fx-font-weight: bold; -fx-font-size: 14px;"); // Larger title font

		// Create internal console if external one not provided
		if (consoleOutput == null) {
			consoleOutput = new TextArea();
		}
		consoleOutput.setEditable(false);
		consoleOutput.setPrefHeight(300); // Set preferred height to make console bigger
		consoleOutput.setMinHeight(200); // Set minimum height
		consoleOutput.setStyle("-fx-font-family: 'Courier New', monospace; -fx-font-size: 12px;"); // Larger, more
																									// readable font
		VBox.setVgrow(consoleOutput, Priority.ALWAYS); // Allow console to grow and fill space

		progressSection.getChildren().addAll(progressLabel, libraryProgressBar, consoleTitle, consoleOutput);

		return progressSection;
	}

	// Library Row class based on MagPepLib structure
	public static class LibraryRow {
		private final SimpleBooleanProperty use;
		private final SimpleStringProperty libName;
		private final SimpleStringProperty pepCount;
		private final SimpleStringProperty protease;
		private final SimpleStringProperty miss;
		private final SimpleStringProperty libType;
		private MagPepLib magPepLib; // Reference to original MagPepLib

		public LibraryRow(boolean use, MagPepLib magPepLib) {
			this.use = new SimpleBooleanProperty(use);
			this.libName = new SimpleStringProperty(magPepLib.getLibName());
			this.pepCount = new SimpleStringProperty(String.valueOf(magPepLib.getPepCount()));
			this.protease = new SimpleStringProperty(magPepLib.getProtease());
			this.miss = new SimpleStringProperty(String.valueOf(magPepLib.getMiss()));
			this.libType = new SimpleStringProperty(magPepLib.isMatchToMag() ? "Microbiome" : "Standalone");
			this.magPepLib = magPepLib;
		}

		public boolean getUse() {
			return use.get();
		}

		public void setUse(boolean value) {
			use.set(value);
		}

		public BooleanProperty useProperty() {
			return use;
		}

		public String getLibName() {
			return libName.get();
		}

		public StringProperty libNameProperty() {
			return libName;
		}

		public String getPepCount() {
			return pepCount.get();
		}

		public StringProperty pepCountProperty() {
			return pepCount;
		}

		public String getProtease() {
			return protease.get();
		}

		public StringProperty proteaseProperty() {
			return protease;
		}

		public String getMiss() {
			return miss.get();
		}

		public StringProperty missProperty() {
			return miss;
		}

		public MagPepLib getMagPepLib() {
			return magPepLib;
		}

		public String getLibType() {
			return libType.get();
		}

		public void setLibType(String value) {
			libType.set(value);
		}

		public StringProperty libTypeProperty() {
			return libType;
		}
	}

	private void browseInputFile() {
		javafx.stage.FileChooser chooser = new javafx.stage.FileChooser();

		RadioButton selected = (RadioButton) inputTypeGroup.getSelectedToggle();
		if (selected.getText().contains("Library")) {
			chooser.setTitle("Select Library File");
			chooser.getExtensionFilters().add(
					new javafx.stage.FileChooser.ExtensionFilter("Library Files", "*.speclib", "*.sptxt", "*.tsv"));
		} else if (selected.getText().contains("Peptide")) {
			chooser.setTitle("Select Peptide List File");
			chooser.getExtensionFilters()
					.add(new javafx.stage.FileChooser.ExtensionFilter("Peptide Files", "*.txt", "*.tsv", "*.csv"));
		} else if (selected.getText().contains("Protein")) {
			chooser.setTitle("Select Protein Database File");
			chooser.getExtensionFilters()
					.add(new javafx.stage.FileChooser.ExtensionFilter("Protein Database Files", "*.fasta", "*.fa"));
		}

		File file = chooser.showOpenDialog(getScene().getWindow());
		if (file != null) {
			inputFileField.setText(file.getAbsolutePath());
		}
	}

	private void generateLibrary() {
		String inputFile = inputFileField.getText();
		String libraryName = libraryNameField.getText();

		if (inputFile == null || inputFile.trim().isEmpty()) {
			showError("Input Required", "Please select an input file.");
			return;
		}

		if (libraryName == null || libraryName.trim().isEmpty()) {
			showError("Name Required", "Please enter a library name.");
			return;
		}

		// Check if MetaParameterDia is available
		if (metaParameter == null) {
			System.out.println("DEBUG: MetaParameterDia is null, attempting to generate from parent panel...");

			// Try to get it from the parent panel
			if (parentPanel != null) {
				System.out.println("DEBUG: Using parent panel to generate MetaParameterDia");
				metaParameter = parentPanel.getCollectedParameters();
				System.out.println("DEBUG: MetaParameterDia generated: " + (metaParameter != null));
			}

			if (metaParameter == null) {
				System.out.println("DEBUG: MetaParameterDia is still null after all attempts");
				showError("Configuration Error",
						"MetaParameterDia is not available. Please ensure all required parameters are set in the main panel.");
				return;
			} else {
				System.out.println("DEBUG: MetaParameterDia successfully obtained");
			}
		}

		if (resourceParameters == null) {
			showError("Configuration Error", "MetaResourceParameters is not set. Please configure the panel properly.");
			return;
		}

		if (currentMagDbItem == null) {
			showError("Database Error", "MAG database item is not set. Please select a database first.");
			return;
		}

		// Determine build type based on selected input type
		RadioButton selected = (RadioButton) inputTypeGroup.getSelectedToggle();
		String sourceType = selected.getText();
		int buildType;

		if (sourceType.contains("Library")) {
			buildType = DiaLibCreateTask.fromLibrary;
		} else if (sourceType.contains("Protein")) {
			buildType = DiaLibCreateTask.fromProtein;
		} else if (sourceType.contains("Peptide")) {
			buildType = DiaLibCreateTask.fromPeptide;
		} else {
			showError("Invalid Input Type", "Please select a valid input type.");
			return;
		}

		// Prepare parameters for DiaLibCreateTask (type can be set per library in the table after creation)
		ArrayList<String> pathList = new ArrayList<>();
		pathList.add(enzymeComboBox.getValue());
		pathList.add(libraryName);
		pathList.add(String.valueOf(missCleavageComboBox.getValue()));
		pathList.add("true"); // default Microbiome; user can change in table

		// Create and run the real DiaLibCreateTask
		DiaLibCreateTask libCreateTask = new DiaLibCreateTask(metaParameter, resourceParameters, libraryProgressBar,
				inputFile, pathList, buildType, LogManager.getLogger(DiaLibraryGenerationPanel.class));

		// Set up task event handlers
		libCreateTask.setOnRunning(e -> {
			Platform.runLater(() -> {
				// Initialize progress tracking
				currentStep = 0;
				updateProgressStep(0, "Starting library generation...");
				libraryProgressBar.setProgress(-1.0);

				if (consoleOutput != null) {
					consoleOutput.appendText("=== Library Generation Started ===\n");
					consoleOutput.appendText("Library name: " + libraryName + "\n");
					consoleOutput.appendText("Input file: " + inputFile + "\n");
					consoleOutput.appendText("Build type: " + sourceType + "\n");
					if (buildType == DiaLibCreateTask.fromProtein) {
						consoleOutput.appendText("Enzyme: " + enzymeComboBox.getValue() + "\n");
						consoleOutput.appendText("Max miss cleavages: " + missCleavageComboBox.getValue() + "\n");
					}
					consoleOutput.appendText("Progress: Step 1/4 - Initializing...\n");
				}
			});
		});

		libCreateTask.setOnSucceeded(e -> {
			Platform.runLater(() -> {
				// Complete progress tracking
				updateProgressStep(totalSteps, "Library generation completed successfully!");
				libraryProgressBar.setProgress(1.0);

				if (consoleOutput != null) {
					consoleOutput.appendText("=== Library Generation Completed ===\n");
					consoleOutput.appendText("✓ Library '" + libraryName + "' generated successfully!\n");
					consoleOutput.appendText("✓ Library has been added to the MAG database.\n");
					consoleOutput.appendText("✓ All steps completed (4/4)\n");
				}

				// Update the configuration file with the new library information
				if (configFile != null && currentMagDbItem != null) {
					try {
						// Determine if this is a MGnify or local MAG database
						boolean isMgnify = currentMagDbItem.getCatalogueID().startsWith("MGYG")
								|| currentMagDbItem.getCatalogueID().contains("mgnify")
								|| currentMagDbItem.getCatalogueID().contains("MGnify");

						// Update the configuration file with the updated MagDbItem
						MagDbConfigIO.updateMagDbItem(currentMagDbItem, isMgnify, configFile);

						if (consoleOutput != null) {
							consoleOutput.appendText("Configuration file updated successfully.\n");
						}

					} catch (Exception ex) {
						if (consoleOutput != null) {
							consoleOutput.appendText(
									"Warning: Failed to update configuration file: " + ex.getMessage() + "\n");
						}
						System.err.println("Failed to update MAG database configuration: " + ex.getMessage());
					}
				} else {
					if (consoleOutput != null) {
						consoleOutput.appendText(
								"Warning: Configuration file or MAG database item not set. Configuration not updated.\n");
					}
				}

				// Refresh the library table
				loadAvailableLibraries();
				notifyParameterUpdate();

				// Clear input fields
				inputFileField.clear();
				libraryNameField.clear();

				// Reset progress after a delay
				Timeline resetTimeline = new Timeline(new KeyFrame(Duration.seconds(3), ev -> {
					libraryProgressBar.setProgress(0);
				}));
				resetTimeline.play();
			});
		});

		libCreateTask.setOnFailed(e -> {
			Platform.runLater(() -> {
				// Reset progress tracking on failure
				updateProgressStep(0, "Library generation failed!");
				libraryProgressBar.setProgress(0);

				Throwable exception = libCreateTask.getException();
				String errorMessage = exception != null ? exception.getMessage() : "Unknown error occurred";

				if (consoleOutput != null) {
					consoleOutput.appendText("=== Library Generation Failed ===\n");
					consoleOutput.appendText("✗ Error at step " + (currentStep + 1) + "/" + totalSteps + "\n");
					consoleOutput.appendText("✗ Error message: " + errorMessage + "\n");
					consoleOutput.appendText("✗ Please check the error details and try again.\n");
				}

				showError("Library Generation Failed", "Failed to generate library at step " + (currentStep + 1) + "/"
						+ totalSteps + ": " + errorMessage);
			});
		});

		libCreateTask.setOnCancelled(e -> {
			Platform.runLater(() -> {
				// Reset progress tracking on cancellation
				updateProgressStep(0, "Library generation was cancelled");
				libraryProgressBar.setProgress(0);

				if (consoleOutput != null) {
					consoleOutput.appendText("=== Library Generation Cancelled ===\n");
					consoleOutput.appendText("⚠ Library generation was cancelled by user.\n");
					consoleOutput.appendText("⚠ Progress: " + currentStep + "/" + totalSteps + " steps completed.\n");
				}
			});
		});

		// Set up progress monitoring from the task
		libCreateTask.messageProperty().addListener((obs, oldMsg, newMsg) -> {
			if (newMsg != null && !newMsg.isEmpty()) {
				Platform.runLater(() -> {
					// Parse progress messages and update step tracking
					if (newMsg.contains("Initializing") || newMsg.contains("Starting")) {
						updateProgressStep(1, "Step 1/4: " + newMsg);
					} else if (newMsg.contains("Processing") || newMsg.contains("Reading")) {
						updateProgressStep(2, "Step 2/4: " + newMsg);
					} else if (newMsg.contains("Building") || newMsg.contains("Creating")) {
						updateProgressStep(3, "Step 3/4: " + newMsg);
					} else if (newMsg.contains("Finalizing") || newMsg.contains("Completing")) {
						updateProgressStep(4, "Step 4/4: " + newMsg);
					}

					if (consoleOutput != null) {
						consoleOutput.appendText("Progress: " + newMsg + "\n");
					}
				});
			}
		});

		// Start the task
		new Thread(libCreateTask).start();
	}

	private void checkAvailableLibraries() {
		if (currentMagDbItem == null) {
			showError("No Database", "Please select a MAG database first.");
			return;
		}

		Task<Void> checkTask = new Task<Void>() {
			@Override
			protected Void call() throws Exception {
				Platform.runLater(() -> {
					if (consoleOutput != null) {
						consoleOutput.appendText("Scanning for available libraries in MAG database folder...\n");
					}
				});

				try {
					// Get the MAG database folder
					File magDbFolder = currentMagDbItem.getCurrentFile();
					if (magDbFolder == null || !magDbFolder.exists()) {
						Platform.runLater(() -> {
							if (consoleOutput != null) {
								consoleOutput.appendText("ERROR: MAG database folder not found: "
										+ (magDbFolder != null ? magDbFolder.getAbsolutePath() : "null") + "\n");
							}
						});
						return null;
					}

					currentMagDbItem.moveLegacyMarker4SpeclibToCanonical();

					// Check the libraries folder
					File librariesFolder = new File(magDbFolder, "libraries");
					if (!librariesFolder.exists()) {
						Platform.runLater(() -> {
							if (consoleOutput != null) {
								consoleOutput.appendText(
										"No libraries folder found at: " + librariesFolder.getAbsolutePath() + "\n");
							}
						});
						return null;
					}

					Platform.runLater(() -> {
						if (consoleOutput != null) {
							consoleOutput.appendText(
									"Scanning libraries folder: " + librariesFolder.getAbsolutePath() + "\n");
						}
					});

					// Scan for library folders
					File[] libraryFolders = librariesFolder.listFiles(File::isDirectory);
					if (libraryFolders == null || libraryFolders.length == 0) {
						Platform.runLater(() -> {
							if (consoleOutput != null) {
								consoleOutput.appendText("No library folders found in libraries directory.\n");
							}
						});
						return null;
					}

					final int[] foundLibraries = { 0 };
					for (File libraryFolder : libraryFolders) {
						String libraryName = libraryFolder.getName();

						File parquetFile = new File(libraryFolder, libraryName + ".parquet");

						if (parquetFile.exists()) {
							Platform.runLater(() -> {
								if (consoleOutput != null) {
									consoleOutput.appendText("Found library: " + libraryName + " (parquet: "
											+ parquetFile.exists() + ")\n");
								}
							});

							// Add the library to the MagDbItem
							System.out.println(
									"DEBUG: Attempting to add library from folder: " + libraryFolder.getAbsolutePath());

							// Check if the required info file exists
							File infoFile = new File(libraryFolder, libraryName + "_info.txt");
							System.out.println("DEBUG: Looking for info file: " + infoFile.getAbsolutePath()
									+ " (exists: " + infoFile.exists() + ")");

							if (infoFile.exists()) {
								currentMagDbItem.addPepLibFromFile(libraryFolder);
								System.out.println("DEBUG: Successfully added library via addPepLibFromFile()");
							} else {
								// Create a basic library entry without the info file
								// We'll use placeholder values since we don't have the info file
								System.out.println("DEBUG: Info file not found, creating basic library entry");
								currentMagDbItem.addPepLib(libraryName, 0, "Trypsin/P", 1);
								System.out.println(
										"DEBUG: Successfully added library via addPepLib() with placeholder values");
							}

							foundLibraries[0]++;
						} else {
							Platform.runLater(() -> {
								if (consoleOutput != null) {
									consoleOutput
											.appendText("Skipping folder (no library files): " + libraryName + "\n");
								}
							});
						}
					}

					currentMagDbItem.registerMarker4SpectralLibraryIfPresent();

					Platform.runLater(() -> {
						if (consoleOutput != null) {
							consoleOutput
									.appendText("Library scan completed. Found " + foundLibraries[0] + " libraries.\n");
						}
						loadAvailableLibraries();
					});

				} catch (Exception e) {
					Platform.runLater(() -> {
						if (consoleOutput != null) {
							consoleOutput.appendText("ERROR during library scan: " + e.getMessage() + "\n");
						}
					});
					System.err.println("Error scanning libraries: " + e.getMessage());
					e.printStackTrace();
				}

				return null;
			}
		};

		new Thread(checkTask).start();
	}

	private void deleteSelectedLibrary() {
		LibraryRow selected = libraryTable.getSelectionModel().getSelectedItem();
		if (selected == null) {
			showError("No Selection", "Please select a library to delete.");
			return;
		}

		// Create confirmation dialog with three options
		Alert alert = new Alert(Alert.AlertType.CONFIRMATION);
		alert.setTitle("Delete Library");
		alert.setHeaderText("Delete library: " + selected.getLibName());
		alert.setContentText("Choose delete option:");

		ButtonType tableOnlyButton = new ButtonType("Remove from table only");
		ButtonType deleteFilesButton = new ButtonType("Delete from table and disk files");
		ButtonType cancelButton = new ButtonType("Cancel", ButtonBar.ButtonData.CANCEL_CLOSE);

		alert.getButtonTypes().setAll(tableOnlyButton, deleteFilesButton, cancelButton);

		Optional<ButtonType> result = alert.showAndWait();
		if (result.isPresent()) {
			if (result.get() == tableOnlyButton) {
				// Remove from table and MAG database only
				if (currentMagDbItem != null) {
					currentMagDbItem.removePepLib(selected.getMagPepLib());
				}
				libraryData.remove(selected);
				notifyParameterUpdate();

				if (consoleOutput != null) {
					consoleOutput.appendText("Library removed from table: " + selected.getLibName() + "\n");
				}
			} else if (result.get() == deleteFilesButton) {
				// Remove from table, MAG database, and delete files
				if (currentMagDbItem != null) {
					currentMagDbItem.removePepLib(selected.getMagPepLib());
				}
				// Delete files from disk
				selected.getMagPepLib().delete();
				libraryData.remove(selected);
				notifyParameterUpdate();

				if (consoleOutput != null) {
					consoleOutput.appendText("Library and files deleted: " + selected.getLibName() + "\n");
				}
			}
			// Cancel - do nothing
		}
	}

	private void loadAvailableLibraries() {
		libraryData.clear();
		if (currentMagDbItem != null) {
			MagPepLib[] magPepLibs = currentMagDbItem.getMagPepLibs();
			System.out.println("DEBUG: Loading libraries from MagDbItem " + currentMagDbItem.getCatalogueID()
					+ " - found " + magPepLibs.length + " libraries");

			for (MagPepLib pepLib : magPepLibs) {
				libraryData.add(new LibraryRow(false, pepLib));
				System.out.println("DEBUG: Added library to table: " + pepLib.getLibName() + " (" + pepLib.getPepCount()
						+ " peptides)");
			}

			if (consoleOutput != null) {
				consoleOutput.appendText("Loaded " + magPepLibs.length + " libraries from MAG database.\n");
			}
		} else {
			System.out.println("DEBUG: Cannot load libraries - currentMagDbItem is null");
			if (consoleOutput != null) {
				consoleOutput.appendText("No MAG database selected - cannot load libraries.\n");
			}
		}
	}

	private Label createSectionTitle(String title) {
		Label label = new Label(title);
		// Use default font instead of "System" to avoid PrismTextLayout issues
		label.setFont(Font.font(null, FontWeight.BOLD, 14));
		label.setTextFill(Color.DARKBLUE);
		return label;
	}

	private void showError(String title, String content) {
		Alert alert = new Alert(Alert.AlertType.ERROR);
		alert.setTitle(title);
		alert.setHeaderText(null);
		alert.setContentText(content);
		alert.showAndWait();
	}

	private void notifyParameterUpdate() {
		if (parameterCallback != null) {
			parameterCallback.onParametersChanged();
		}
	}

	// Public methods for accessing selected libraries
	public MagPepLib[] getSelectedLibraries() {
		return libraryData.stream().filter(LibraryRow::getUse).map(LibraryRow::getMagPepLib).toArray(MagPepLib[]::new);
	}

	public String[] getSelectedLibraryNames() {
		System.out.println("DEBUG: getSelectedLibraryNames() called");
		System.out.println("DEBUG: Total libraries in table: " + libraryData.size());

		String[] selectedLibs = libraryData.stream().filter(LibraryRow::getUse).map(LibraryRow::getLibName)
				.toArray(String[]::new);

		System.out.println("DEBUG: Selected libraries count: " + selectedLibs.length);
		for (int i = 0; i < selectedLibs.length; i++) {
			System.out.println("DEBUG: Selected library " + i + ": " + selectedLibs[i]);
		}

		// Also debug all libraries and their selection status
		for (int i = 0; i < libraryData.size(); i++) {
			LibraryRow row = libraryData.get(i);
			System.out.println("DEBUG: Library " + i + ": " + row.getLibName() + " (selected: " + row.getUse() + ")");
		}

		return selectedLibs;
	}

	/**
	 * Public method to refresh the library table from the MagDbItem This can be
	 * called externally to update the table when libraries are added/removed
	 */
	public void refreshLibraryTable() {
		System.out.println("DEBUG: Refreshing library table");
		loadAvailableLibraries();
	}

	/**
	 * Test method to manually select the first library for debugging
	 */
	public void testSelectFirstLibrary() {
		if (libraryData.size() > 0) {
			LibraryRow firstRow = libraryData.get(0);
			firstRow.setUse(true);
			System.out.println("DEBUG: Manually selected first library: " + firstRow.getLibName());
			notifyParameterUpdate();
		} else {
			System.out.println("DEBUG: No libraries available to select");
		}
	}

	/**
	 * Update progress step with clear indication
	 */
	private void updateProgressStep(int step, String message) {
		currentStep = step;

		if (progressLabel != null) {
			if (step == 0) {
				progressLabel.setText("Ready to generate library");
				progressLabel.setStyle("-fx-font-weight: bold; -fx-font-size: 12px; -fx-text-fill: #333333;");
			} else if (step <= totalSteps) {
				progressLabel.setText("Step " + step + "/" + totalSteps + ": " + message);
				progressLabel.setStyle("-fx-font-weight: bold; -fx-font-size: 12px; -fx-text-fill: #0066CC;");
			} else {
				progressLabel.setText("✓ " + message);
				progressLabel.setStyle("-fx-font-weight: bold; -fx-font-size: 12px; -fx-text-fill: #009900;");
			}
		}

		System.out.println("DEBUG: Progress updated - Step " + step + "/" + totalSteps + ": " + message);
	}
}
