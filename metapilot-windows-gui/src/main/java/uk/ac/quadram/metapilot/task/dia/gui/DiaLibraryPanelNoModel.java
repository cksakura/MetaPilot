package uk.ac.quadram.metapilot.task.dia.gui;

import javafx.application.Platform;
import javafx.geometry.Insets;
import javafx.geometry.Pos;
import javafx.scene.control.*;
import uk.ac.quadram.metapilot.task.common.gui.ParameterSection;
import javafx.scene.layout.*;
import javafx.collections.FXCollections;
import javafx.collections.ObservableList;
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
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.Set;
import java.util.concurrent.CountDownLatch;
import java.util.concurrent.atomic.AtomicBoolean;

import org.apache.logging.log4j.LogManager;

import javafx.animation.KeyFrame;
import javafx.animation.Timeline;
import javafx.util.Duration;
import uk.ac.quadram.metapilot.task.dia.DiaLibCreateTask;
import uk.ac.quadram.metapilot.task.dia.DiaMarkerCoreTask;
import uk.ac.quadram.metapilot.task.dia.par.MetaParameterDia;
import uk.ac.quadram.metapilot.resources.MetaResourceParameters;
import uk.ac.quadram.metapilot.task.mag.MagDbConfigIO;

/**
 * Library panel for DIA Normal workflow (no ML model).
 * Split layout: controls section and table section can be placed in different columns
 * for better view when there are many library items.
 */
public class DiaLibraryPanelNoModel extends VBox {

	public interface ParameterUpdateCallback {
		void onParametersChanged();
	}

	public interface ParentPanel {
		MetaParameterDia collectParametersFromUI();
	}

	private TableView<LibraryRow> libraryTable;
	private ObservableList<LibraryRow> libraryData = FXCollections.observableArrayList();

	private ToggleGroup inputTypeGroup;
	private TextField inputFileField;
	private ComboBox<String> enzymeComboBox;
	private ComboBox<Integer> missCleavageComboBox;
	private TextField libraryNameField;
	private VBox enzymeBox;

	private MagDbItem currentMagDbItem;
	private TextArea consoleOutput;
	private Button generateLibraryButton;
	private ProgressBar libraryProgressBar;
	private Label progressLabel;

	private int currentStep = 0;
	private int totalSteps = 4;

	private ParameterUpdateCallback parameterCallback;
	private ParentPanel parentPanel;

	private MetaParameterDia metaParameter;
	private MetaResourceParameters resourceParameters;
	private File configFile;

	private TitledPane controlsSection;
	private TitledPane tableSection;

	public DiaLibraryPanelNoModel(MetaParameterDia metaParameter, MetaResourceParameters resourceParameters,
			File configFile) {
		this.metaParameter = metaParameter;
		this.resourceParameters = resourceParameters;
		this.configFile = configFile;
		if (metaParameter != null && metaParameter.getMagDbItem() != null) {
			this.currentMagDbItem = metaParameter.getMagDbItem();
		}
		initializeUI();
		loadAvailableLibraries();
	}

	public void setParameterUpdateCallback(ParameterUpdateCallback callback) {
		this.parameterCallback = callback;
	}

	public void setParentPanel(ParentPanel parentPanel) {
		this.parentPanel = parentPanel;
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
		applyLibrarySelectionFromParameters();
	}

	public void setResourceParameters(MetaResourceParameters resourceParameters) {
		this.resourceParameters = resourceParameters;
	}

	public void setConfigFile(File configFile) {
		this.configFile = configFile;
	}

	public void updateMetaParameter(MetaParameterDia metaParameter) {
		this.metaParameter = metaParameter;
		applyLibrarySelectionFromParameters();
	}

	/** Returns the controls section (input type, file, enzyme, generate button, progress) for left column. */
	public TitledPane getControlsSection() {
		return controlsSection;
	}

	/** Returns the library table section for right column - better view for many items. */
	public TitledPane getTableSection() {
		return tableSection;
	}

	public MagPepLib[] getSelectedLibraries() {
		return libraryData.stream().filter(LibraryRow::getUse).map(LibraryRow::getMagPepLib).toArray(MagPepLib[]::new);
	}

	public String[] getSelectedLibraryNames() {
		return libraryData.stream().filter(LibraryRow::getUse).map(LibraryRow::getLibName).toArray(String[]::new);
	}

	public void refreshLibraryTable() {
		loadAvailableLibraries();
	}

	private void initializeUI() {
		controlsSection = createControlsSection();
		tableSection = createTableSection();
	}

	private TitledPane createControlsSection() {
		VBox section = new VBox(10);


		section.setMinWidth(350);


		HBox inputTypeBox = createInputTypeSelection();
		HBox fileBox = createFileSelection();
		enzymeBox = createEnzymeParameters();
		HBox generationBox = createGenerationControls();
		VBox progressBox = createProgressSection();

		section.getChildren().addAll(inputTypeBox, fileBox, enzymeBox, generationBox, progressBox);
		VBox.setVgrow(progressBox, Priority.ALWAYS);
		return ParameterSection.create("Spectral Library Generation", section);
	}

	private TitledPane createTableSection() {
		VBox section = new VBox(10);
		section.setMinHeight(300);
		VBox.setVgrow(section, Priority.ALWAYS);

		HBox headerBox = new HBox(10);
		headerBox.setAlignment(Pos.CENTER_LEFT);
		Label availableLabel = new Label("Select libraries for identification");
		availableLabel.setStyle("-fx-font-weight: bold;");

		Button checkAvailableButton = new Button("Check Available");
		checkAvailableButton.setTooltip(new Tooltip("Auto-detect libraries in MAG database folder"));
		checkAvailableButton.setOnAction(e -> checkAvailableLibraries());

		Button deleteButton = new Button("Delete Selected");
		deleteButton.setOnAction(e -> deleteSelectedLibrary());

		headerBox.getChildren().addAll(availableLabel, checkAvailableButton, deleteButton);

		libraryTable = new TableView<>(libraryData);
		libraryTable.setMinHeight(150);
		libraryTable.setPrefHeight(350);
		libraryTable.setMaxHeight(Double.MAX_VALUE);
		libraryTable.setEditable(true);
		libraryTable.setColumnResizePolicy(TableView.UNCONSTRAINED_RESIZE_POLICY);
		VBox.setVgrow(libraryTable, Priority.ALWAYS);

		TableColumn<LibraryRow, Boolean> useCol = new TableColumn<>("Use");
		useCol.setCellValueFactory(new PropertyValueFactory<>("use"));
		useCol.setCellFactory(column -> new TableCell<LibraryRow, Boolean>() {
			private final CheckBox checkBox = new CheckBox();
			{
				checkBox.setOnAction(e -> {
					LibraryRow row = getTableRow().getItem();
					if (row != null) {
						row.setUse(checkBox.isSelected());
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
		});
		useCol.setMinWidth(40);
		useCol.setPrefWidth(50);

		TableColumn<LibraryRow, String> nameCol = new TableColumn<>("Library Name");
		nameCol.setCellValueFactory(new PropertyValueFactory<>("libName"));
		nameCol.setMinWidth(120);

		TableColumn<LibraryRow, String> pepCountCol = new TableColumn<>("Peptide Count");
		pepCountCol.setCellValueFactory(new PropertyValueFactory<>("pepCount"));
		pepCountCol.setMinWidth(80);

		TableColumn<LibraryRow, String> proteaseCol = new TableColumn<>("Protease");
		proteaseCol.setCellValueFactory(new PropertyValueFactory<>("protease"));
		proteaseCol.setMinWidth(80);

		TableColumn<LibraryRow, String> missCol = new TableColumn<>("Miss Cleavages");
		missCol.setCellValueFactory(new PropertyValueFactory<>("miss"));
		missCol.setMinWidth(80);

		TableColumn<LibraryRow, String> typeCol = new TableColumn<>("Type");
		typeCol.setCellValueFactory(new PropertyValueFactory<>("libType"));
		typeCol.setMinWidth(100);
		typeCol.setPrefWidth(140);
		typeCol.setCellFactory(column -> new TableCell<LibraryRow, String>() {
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
		});

		libraryTable.getColumns().addAll(useCol, nameCol, pepCountCol, proteaseCol, missCol, typeCol);

		section.getChildren().addAll(headerBox, libraryTable);
		return ParameterSection.create("Available Libraries", section);
	}

	private HBox createInputTypeSelection() {
		HBox box = new HBox(8);
		Label inputLabel = new Label("Input Type:");
		inputLabel.setStyle("-fx-font-weight: bold;");

		inputTypeGroup = new ToggleGroup();
		RadioButton libraryRadio = new RadioButton("Existing Library File");
		RadioButton peptideListRadio = new RadioButton("Peptide List File");
		RadioButton proteinDbRadio = new RadioButton("Protein Database File");

		libraryRadio.setToggleGroup(inputTypeGroup);
		peptideListRadio.setToggleGroup(inputTypeGroup);
		proteinDbRadio.setToggleGroup(inputTypeGroup);
		libraryRadio.setSelected(true);

		inputTypeGroup.selectedToggleProperty().addListener((obs, oldToggle, newToggle) -> {
			boolean enableEnzyme = newToggle == proteinDbRadio;
			enzymeComboBox.setDisable(!enableEnzyme);
			missCleavageComboBox.setDisable(!enableEnzyme);
		});

		box.getChildren().addAll(inputLabel, libraryRadio, peptideListRadio, proteinDbRadio);
		return box;
	}

	private HBox createFileSelection() {
		HBox box = new HBox(10);
		box.setAlignment(Pos.CENTER_LEFT);
		Label fileLabel = new Label("Input File:");
		fileLabel.setStyle("-fx-font-weight: bold;");
		inputFileField = new TextField();
		inputFileField.setPromptText("Select input file");
		inputFileField.setPrefWidth(250);
		inputFileField.setEditable(false);
		Button browseFileButton = new Button("Browse");
		browseFileButton.setOnAction(e -> browseInputFile());
		box.getChildren().addAll(fileLabel, inputFileField, browseFileButton);
		return box;
	}

	private VBox createEnzymeParameters() {
		VBox box = new VBox(8);
		box.setStyle("-fx-background-color: #f0f0f0; -fx-padding: 10; -fx-background-radius: 3;");
		Label enzymeTitle = new Label("Enzyme Parameters (for Protein Database):");
		enzymeTitle.setStyle("-fx-font-weight: bold;");

		GridPane enzymeGrid = new GridPane();
		enzymeGrid.setHgap(10);
		enzymeGrid.setVgap(5);

		Label enzymeNameLabel = new Label("Enzyme:");
		enzymeComboBox = new ComboBox<>(
				FXCollections.observableArrayList("Trypsin/P", "Trypsin", "Lys-C", "Chymotrypsin", "AspN", "GluC"));
		enzymeComboBox.setValue("Trypsin/P");
		enzymeComboBox.setPrefWidth(120);

		Label missCleavageLabel = new Label("Max Miss Cleavages:");
		missCleavageComboBox = new ComboBox<>();
		missCleavageComboBox.getItems().addAll(0, 1, 2);
		missCleavageComboBox.setValue(1);
		missCleavageComboBox.setPrefWidth(80);

		enzymeGrid.add(enzymeNameLabel, 0, 0);
		enzymeGrid.add(enzymeComboBox, 1, 0);
		enzymeGrid.add(missCleavageLabel, 0, 1);
		enzymeGrid.add(missCleavageComboBox, 1, 1);

		box.getChildren().addAll(enzymeTitle, enzymeGrid);
		enzymeComboBox.setDisable(true);
		missCleavageComboBox.setDisable(true);

		return box;
	}

	private HBox createGenerationControls() {
		HBox box = new HBox(10);
		box.setAlignment(Pos.CENTER_LEFT);
		libraryNameField = new TextField();
		libraryNameField.setPromptText("Library name");
		libraryNameField.setPrefWidth(150);

		if (generateLibraryButton == null) {
			generateLibraryButton = new Button("Generate Library");
			generateLibraryButton.setPrefWidth(130);
			generateLibraryButton.setOnAction(e -> generateLibrary());
		}

		box.getChildren().addAll(new Label("Library Name:"), libraryNameField, generateLibraryButton);
		return box;
	}

	private VBox createProgressSection() {
		VBox box = new VBox(8);
		box.setPadding(new Insets(10));
		box.setStyle("-fx-background-color: #f8f8f8; -fx-background-radius: 5;");

		libraryProgressBar = new ProgressBar(0);
		libraryProgressBar.setMaxWidth(Double.MAX_VALUE);
		libraryProgressBar.setPrefHeight(25);
		libraryProgressBar.setStyle("-fx-accent: #4CAF50;");

		progressLabel = new Label("Ready to generate library");
		progressLabel.setStyle("-fx-font-weight: bold; -fx-font-size: 12px; -fx-text-fill: #333333;");
		progressLabel.setMaxWidth(Double.MAX_VALUE);

		if (consoleOutput == null) {
			consoleOutput = new TextArea();
		}
		consoleOutput.setEditable(false);
		consoleOutput.setMinHeight(120);
		consoleOutput.setPrefHeight(200);
		consoleOutput.setMaxHeight(Double.MAX_VALUE);
		consoleOutput.setStyle("-fx-font-family: 'Courier New', monospace; -fx-font-size: 11px;");
		Label consoleTitle = new Label("Console Output:");
		consoleTitle.setStyle("-fx-font-weight: bold; -fx-font-size: 12px;");

		box.getChildren().addAll(progressLabel, libraryProgressBar, consoleTitle, consoleOutput);
		VBox.setVgrow(consoleOutput, Priority.ALWAYS);
		return box;
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

		File file = chooser.showOpenDialog(getScene() != null ? getScene().getWindow() : null);
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

		if (metaParameter == null && parentPanel != null) {
			metaParameter = parentPanel.collectParametersFromUI();
		}
		if (metaParameter == null) {
			showError("Configuration Error", "Parameters not available. Please ensure all required parameters are set.");
			return;
		}
		if (resourceParameters == null) {
			showError("Configuration Error", "MetaResourceParameters is not set.");
			return;
		}
		if (currentMagDbItem == null) {
			showError("Database Error", "MAG database item is not set. Please select a database first.");
			return;
		}

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

		ArrayList<String> pathList = new ArrayList<>();
		pathList.add(enzymeComboBox.getValue());
		pathList.add(libraryName);
		pathList.add(String.valueOf(missCleavageComboBox.getValue()));
		pathList.add("true"); // default Microbiome; user can change in table

		DiaLibCreateTask libCreateTask = new DiaLibCreateTask(metaParameter, resourceParameters, libraryProgressBar,
				inputFile, pathList, buildType, LogManager.getLogger(DiaLibraryPanelNoModel.class));

		libCreateTask.setOnRunning(e -> Platform.runLater(() -> {
			currentStep = 0;
			if (generateLibraryButton != null) {
				generateLibraryButton.setDisable(true);
			}
			if (libraryProgressBar != null) {
				libraryProgressBar.setProgress(-1.0); // indeterminate while tasks initialize
			}
			updateProgressStep(0, "Starting library generation...");
			if (consoleOutput != null) {
				consoleOutput.appendText("=== Library Generation Started ===\n");
				consoleOutput.appendText("Library name: " + libraryName + "\n");
				consoleOutput.appendText("Input file: " + inputFile + "\n");
			}
		}));

		libCreateTask.setOnSucceeded(e -> Platform.runLater(() -> {
			Boolean taskResult = libCreateTask.getValue();
			if (taskResult == null || !taskResult) {
				updateProgressStep(0, "Library generation failed!");
				if (libraryProgressBar != null) {
					libraryProgressBar.setProgress(0);
				}
				if (generateLibraryButton != null) {
					generateLibraryButton.setDisable(false);
				}
				String taskFailureReason = libCreateTask.getLastFailureReason();
				if (taskFailureReason == null || taskFailureReason.trim().isEmpty()) {
					taskFailureReason = "DiaLibCreateTask returned false without an explicit reason.";
				}
				if (consoleOutput != null) {
					consoleOutput.appendText("=== Library Generation Failed ===\n");
					consoleOutput.appendText("✗ DiaLibCreateTask failed: " + taskFailureReason + "\n");
				}
				showError("Library Generation Failed", "DiaLibCreateTask failed:\n" + taskFailureReason);
				return;
			}

			String outputValidationError = validateGeneratedLibraryOutput(libraryName);
			if (outputValidationError != null) {
				updateProgressStep(0, "Library generation failed (output validation)!");
				if (libraryProgressBar != null) {
					libraryProgressBar.setProgress(0);
				}
				if (generateLibraryButton != null) {
					generateLibraryButton.setDisable(false);
				}
				if (consoleOutput != null) {
					consoleOutput.appendText("=== Library Generation Failed ===\n");
					consoleOutput.appendText("✗ Library output validation failed.\n");
					consoleOutput.appendText(outputValidationError);
				}
				showError("Library Generation Failed",
						"Library task completed but output files were not generated correctly.\n\n" + outputValidationError);
				return;
			}

			updateProgressStep(totalSteps, "Library generation completed successfully!");
			if (libraryProgressBar != null) {
				libraryProgressBar.setProgress(1.0);
			}
			if (generateLibraryButton != null) {
				generateLibraryButton.setDisable(false);
			}
			if (consoleOutput != null) {
				consoleOutput.appendText("=== Library Generation Completed ===\n");
				consoleOutput.appendText("✓ Library '" + libraryName + "' generated successfully!\n");
			}
			if (configFile != null && currentMagDbItem != null) {
				try {
					boolean isMgnify = currentMagDbItem.getCatalogueID().startsWith("MGYG")
							|| currentMagDbItem.getCatalogueID().contains("mgnify")
							|| currentMagDbItem.getCatalogueID().contains("MGnify");
					MagDbConfigIO.updateMagDbItem(currentMagDbItem, isMgnify, configFile);
				} catch (Exception ex) {
					if (consoleOutput != null) {
						consoleOutput.appendText("Warning: Failed to update configuration: " + ex.getMessage() + "\n");
					}
				}
			}
			loadAvailableLibraries();
			notifyParameterUpdate();
			inputFileField.clear();
			libraryNameField.clear();
			Timeline resetTimeline = new Timeline(new KeyFrame(Duration.seconds(3), ev -> libraryProgressBar.setProgress(0)));
			resetTimeline.play();
		}));

		libCreateTask.setOnFailed(e -> Platform.runLater(() -> {
			updateProgressStep(0, "Library generation failed!");
			if (libraryProgressBar != null) {
				libraryProgressBar.setProgress(0);
			}
			if (generateLibraryButton != null) {
				generateLibraryButton.setDisable(false);
			}
			Throwable exception = libCreateTask.getException();
			String errorMessage = exception != null ? exception.getMessage() : "Unknown error occurred";
			if (consoleOutput != null) {
				consoleOutput.appendText("=== Library Generation Failed ===\n");
				consoleOutput.appendText("✗ Error: " + errorMessage + "\n");
			}
			showError("Library Generation Failed", "Failed to generate library: " + errorMessage);
		}));

		libCreateTask.setOnCancelled(e -> Platform.runLater(() -> {
			updateProgressStep(0, "Library generation was cancelled");
			if (libraryProgressBar != null) {
				libraryProgressBar.setProgress(0);
			}
			if (generateLibraryButton != null) {
				generateLibraryButton.setDisable(false);
			}
		}));

		libCreateTask.messageProperty().addListener((obs, oldMsg, newMsg) -> {
			if (newMsg != null && !newMsg.isEmpty()) {
				Platform.runLater(() -> {
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

		new Thread(libCreateTask).start();
	}

	/**
	 * Verifies expected generated library output exists; returns null if valid, otherwise a detailed reason.
	 */
	private String validateGeneratedLibraryOutput(String libraryName) {
		if (currentMagDbItem == null) {
			return "Reason: current MAG database is null.\n";
		}
		File magDbFolder = currentMagDbItem.getCurrentFile();
		if (magDbFolder == null) {
			return "Reason: MAG database folder is null.\n";
		}
		if (!magDbFolder.exists()) {
			return "Reason: MAG database folder does not exist: " + magDbFolder.getAbsolutePath() + "\n";
		}

		File librariesFolder = new File(magDbFolder, "libraries");
		if (!librariesFolder.exists()) {
			return "Reason: libraries folder does not exist: " + librariesFolder.getAbsolutePath() + "\n";
		}
		if (!librariesFolder.isDirectory()) {
			return "Reason: libraries path is not a directory: " + librariesFolder.getAbsolutePath() + "\n";
		}

		File libraryFolder = new File(librariesFolder, libraryName);
		if (!libraryFolder.exists()) {
			return "Reason: generated library folder was not created: " + libraryFolder.getAbsolutePath() + "\n";
		}
		if (!libraryFolder.isDirectory()) {
			return "Reason: generated library path is not a directory: " + libraryFolder.getAbsolutePath() + "\n";
		}

		File[] contents = libraryFolder.listFiles();
		if (contents == null) {
			return "Reason: cannot read generated library folder: " + libraryFolder.getAbsolutePath() + "\n";
		}
		if (contents.length == 0) {
			return "Reason: generated library folder is empty: " + libraryFolder.getAbsolutePath() + "\n";
		}

		File parquetFile = new File(libraryFolder, libraryName + ".parquet");
		if (!parquetFile.exists()) {
			StringBuilder sb = new StringBuilder();
			sb.append("Reason: expected parquet library file was not found in: ").append(libraryFolder.getAbsolutePath()).append("\n");
			sb.append("Expected: ").append(parquetFile.getName()).append("\n");
			sb.append("Actual files: ");
			for (int i = 0; i < contents.length; i++) {
				if (i > 0) {
					sb.append(", ");
				}
				sb.append(contents[i].getName());
			}
			sb.append("\n");
			return sb.toString();
		}
		return null;
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
						consoleOutput.appendText("Scanning for available libraries...\n");
					}
				});

				try {
					File magDbFolder = currentMagDbItem.getCurrentFile();
					if (magDbFolder == null || !magDbFolder.exists()) {
						Platform.runLater(() -> {
							if (consoleOutput != null) {
								consoleOutput.appendText("ERROR: MAG database folder not found.\n");
								consoleOutput.appendText("MAG database path: "
										+ (magDbFolder == null ? "<null>" : magDbFolder.getAbsolutePath()) + "\n");
							}
						});
						return null;
					}

					currentMagDbItem.moveLegacyMarker4SpeclibToCanonical();

					File marker4Speclib = currentMagDbItem.getMarker4LibFile();
					boolean marker4Missing = marker4Speclib == null || !marker4Speclib.exists();
					if (marker4Missing && resourceParameters != null) {
						CountDownLatch latch = new CountDownLatch(1);
						AtomicBoolean userWantsGeneration = new AtomicBoolean(false);
						Platform.runLater(() -> {
							try {
								Alert confirm = new Alert(Alert.AlertType.CONFIRMATION);
								confirm.setTitle("Marker4 spectral library");
								confirm.setHeaderText("No Marker4 library found");
								confirm.setContentText(
										"Generate the Marker4 spectral library now from your MAG database FASTA?\n\n"
												+ "This uses DIA-NN (predicted spectra) and may take several minutes. "
												+ "You can skip and run a Marker DIA task later instead.");
								confirm.getButtonTypes().setAll(ButtonType.YES, ButtonType.NO);
								Optional<ButtonType> answer = confirm.showAndWait();
								userWantsGeneration.set(answer.isPresent() && answer.get() == ButtonType.YES);
							} finally {
								latch.countDown();
							}
						});
						latch.await();

						if (userWantsGeneration.get()) {
							MetaParameterDia diaPar = metaParameter;
							if (parentPanel != null) {
								try {
									MetaParameterDia collected = parentPanel.collectParametersFromUI();
									if (collected != null) {
										diaPar = collected;
									}
								} catch (Exception ignore) {
									// keep metaParameter
								}
							}
							if (diaPar == null) {
								Platform.runLater(
										() -> showError("Parameters", "Cannot generate Marker4: workflow parameters are not available."));
							} else {
								Platform.runLater(() -> {
									if (consoleOutput != null) {
										consoleOutput.appendText("Generating Marker4 spectral library...\n");
									}
								});
								boolean genOk = DiaMarkerCoreTask.generateMarker4LibraryOnly(diaPar, resourceParameters,
										currentMagDbItem, libraryProgressBar, null);
								if (!genOk) {
									Platform.runLater(() -> {
										if (consoleOutput != null) {
											consoleOutput.appendText(
													"ERROR: Marker4 library generation failed. Check that DIA-NN is configured, "
															+ "Marker4 FASTA exists for this MAG database, and see the log.\n");
										}
										showError("Marker4 generation failed",
												"Could not build the Marker4 library. Check DIA-NN in resource settings, "
														+ "ensure the MAG database has Marker4 proteins, and review the console log.");
									});
								} else {
									Platform.runLater(() -> {
										if (consoleOutput != null) {
											consoleOutput.appendText("Marker4 library generated successfully.\n");
										}
									});
								}
							}
						}
					}

					File librariesFolder = new File(magDbFolder, MagDbItem.LIBRARIES_FOLDER_NAME);
					Platform.runLater(() -> {
						if (consoleOutput != null) {
							consoleOutput.appendText("MAG database folder: " + magDbFolder.getAbsolutePath() + "\n");
							consoleOutput.appendText("Expected libraries folder: " + librariesFolder.getAbsolutePath() + "\n");
						}
					});
					if (!librariesFolder.exists()) {
						Platform.runLater(() -> {
							if (consoleOutput != null) {
								consoleOutput.appendText("No libraries folder found.\n");
								consoleOutput.appendText("Reason details: exists=" + librariesFolder.exists()
										+ ", isDirectory=" + librariesFolder.isDirectory() + ", canRead="
										+ librariesFolder.canRead() + "\n");
							}
						});
						return null;
					}

					File[] libraryFolders = librariesFolder.listFiles(File::isDirectory);
					if (libraryFolders == null || libraryFolders.length == 0) {
						Platform.runLater(() -> {
							if (consoleOutput != null) {
								consoleOutput.appendText("No library folders found in: " + librariesFolder.getAbsolutePath() + "\n");
								consoleOutput.appendText("Reason details: exists=" + librariesFolder.exists()
										+ ", isDirectory=" + librariesFolder.isDirectory() + ", canRead="
										+ librariesFolder.canRead() + ", listFilesNull=" + (libraryFolders == null) + "\n");
							}
						});
						return null;
					}

					int[] foundLibraries = { 0 };
					for (File libraryFolder : libraryFolders) {
						String libraryName = libraryFolder.getName();
						File parquetFile = new File(libraryFolder, libraryName + ".parquet");

						if (parquetFile.exists()) {
							File infoFile = new File(libraryFolder, libraryName + "_info.txt");
							if (infoFile.exists()) {
								currentMagDbItem.addPepLibFromFile(libraryFolder);
							} else {
								currentMagDbItem.addPepLib(libraryName, 0, "Trypsin/P", 1);
							}
							foundLibraries[0]++;
						}
					}

					currentMagDbItem.registerMarker4SpectralLibraryIfPresent();

					final int count = foundLibraries[0];
					Platform.runLater(() -> {
						if (consoleOutput != null) {
							consoleOutput.appendText("Library scan completed. Found " + count + " libraries.\n");
						}
						loadAvailableLibraries();
					});

				} catch (Exception e) {
					Platform.runLater(() -> {
						if (consoleOutput != null) {
							consoleOutput.appendText("ERROR during library scan: " + e.getMessage() + "\n");
						}
					});
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
				if (currentMagDbItem != null) {
					currentMagDbItem.removePepLib(selected.getMagPepLib());
				}
				libraryData.remove(selected);
				notifyParameterUpdate();
			} else if (result.get() == deleteFilesButton) {
				if (currentMagDbItem != null) {
					currentMagDbItem.removePepLib(selected.getMagPepLib());
				}
				selected.getMagPepLib().delete();
				libraryData.remove(selected);
				notifyParameterUpdate();
			}
		}
	}

	private void loadAvailableLibraries() {
		libraryData.clear();
		if (currentMagDbItem != null) {
			MagPepLib[] magPepLibs = currentMagDbItem.getMagPepLibs();
			for (MagPepLib pepLib : magPepLibs) {
				libraryData.add(new LibraryRow(false, pepLib));
			}
			applyLibrarySelectionFromParameters();
			if (consoleOutput != null) {
				consoleOutput.appendText("Loaded " + magPepLibs.length + " libraries from MAG database.\n");
			}
		}
	}

	/**
	 * Apply saved library names (MetaParameterDia.library) to table checkboxes.
	 * This is critical for round-trip persistence: file -> UI -> file.
	 */
	private void applyLibrarySelectionFromParameters() {
		if (metaParameter == null || libraryData == null || libraryData.isEmpty()) {
			return;
		}
		String[] selected = metaParameter.getLibrary();
		Set<String> selectedNames = new HashSet<>();
		if (selected != null && selected.length > 0) {
			selectedNames.addAll(Arrays.asList(selected));
		}
		for (LibraryRow row : libraryData) {
			row.setUse(selectedNames.contains(row.getLibName()));
		}
		if (libraryTable != null) {
			libraryTable.refresh();
		}
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
	}

	public static class LibraryRow {
		private final SimpleBooleanProperty use;
		private final SimpleStringProperty libName;
		private final SimpleStringProperty pepCount;
		private final SimpleStringProperty protease;
		private final SimpleStringProperty miss;
		private final SimpleStringProperty libType;
		private MagPepLib magPepLib;

		public LibraryRow(boolean use, MagPepLib magPepLib) {
			this.use = new SimpleBooleanProperty(use);
			this.libName = new SimpleStringProperty(magPepLib.getLibName());
			this.pepCount = new SimpleStringProperty(String.valueOf(magPepLib.getPepCount()));
			this.protease = new SimpleStringProperty(magPepLib.getProtease());
			this.miss = new SimpleStringProperty(String.valueOf(magPepLib.getMiss()));
			this.libType = new SimpleStringProperty(magPepLib.isMatchToMag() ? "Microbiome" : "Standalone");
			this.magPepLib = magPepLib;
		}

		public boolean getUse() { return use.get(); }
		public void setUse(boolean value) { use.set(value); }
		public BooleanProperty useProperty() { return use; }
		public String getLibName() { return libName.get(); }
		public StringProperty libNameProperty() { return libName; }
		public String getPepCount() { return pepCount.get(); }
		public StringProperty pepCountProperty() { return pepCount; }
		public String getProtease() { return protease.get(); }
		public StringProperty proteaseProperty() { return protease; }
		public String getMiss() { return miss.get(); }
		public StringProperty missProperty() { return miss; }
		public String getLibType() { return libType.get(); }
		public void setLibType(String value) { libType.set(value); }
		public StringProperty libTypeProperty() { return libType; }
		public MagPepLib getMagPepLib() { return magPepLib; }
	}
}
