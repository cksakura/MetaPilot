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
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.Set;

import org.apache.logging.log4j.LogManager;

import javafx.animation.KeyFrame;
import javafx.animation.Timeline;
import javafx.util.Duration;
import uk.ac.quadram.metapilot.task.dia.DiaLibCreateTask;
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

	private VBox controlsSection;
	private VBox tableSection;

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
	public VBox getControlsSection() {
		return controlsSection;
	}

	/** Returns the library table section for right column - better view for many items. */
	public VBox getTableSection() {
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

	private VBox createControlsSection() {
		VBox section = new VBox(10);
		section.setPadding(new Insets(10));
		section.setStyle("-fx-background-color: #fff8dc; -fx-background-radius: 5;");
		section.setMinWidth(350);

		Label titleLabel = createSectionTitle("Spectral Library Generation");

		HBox inputTypeBox = createInputTypeSelection();
		HBox fileBox = createFileSelection();
		enzymeBox = createEnzymeParameters();
		HBox generationBox = createGenerationControls();
		VBox progressBox = createProgressSection();

		section.getChildren().addAll(titleLabel, inputTypeBox, fileBox, enzymeBox, generationBox, progressBox);
		VBox.setVgrow(progressBox, Priority.ALWAYS);
		return section;
	}

	private VBox createTableSection() {
		VBox section = new VBox(10);
		section.setPadding(new Insets(10));
		section.setStyle("-fx-background-color: #f8f8f8; -fx-background-radius: 5;");
		section.setMinHeight(300);
		VBox.setVgrow(section, Priority.ALWAYS);

		HBox headerBox = new HBox(10);
		headerBox.setAlignment(Pos.CENTER_LEFT);
		Label availableLabel = new Label("Available Libraries (select for identification):");
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
		return section;
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
							}
						});
						return null;
					}

					File librariesFolder = new File(magDbFolder, "libraries");
					if (!librariesFolder.exists()) {
						Platform.runLater(() -> {
							if (consoleOutput != null) {
								consoleOutput.appendText("No libraries folder found.\n");
							}
						});
						return null;
					}

					File[] libraryFolders = librariesFolder.listFiles(File::isDirectory);
					if (libraryFolders == null || libraryFolders.length == 0) {
						Platform.runLater(() -> {
							if (consoleOutput != null) {
								consoleOutput.appendText("No library folders found.\n");
							}
						});
						return null;
					}

					int[] foundLibraries = { 0 };
					for (File libraryFolder : libraryFolders) {
						String libraryName = libraryFolder.getName();
						File speclibFile = new File(libraryFolder, libraryName + ".speclib");
						File tsvFile = new File(libraryFolder, libraryName + ".tsv");
						File parquetFile = new File(libraryFolder, libraryName + ".parquet");

						if (speclibFile.exists() || tsvFile.exists() || parquetFile.exists()) {
							File infoFile = new File(libraryFolder, libraryName + "_info.txt");
							if (infoFile.exists()) {
								currentMagDbItem.addPepLibFromFile(libraryFolder);
							} else {
								currentMagDbItem.addPepLib(libraryName, 0, "Trypsin/P", 1);
							}
							foundLibraries[0]++;
						}
					}

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

	private Label createSectionTitle(String title) {
		Label label = new Label(title);
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

	private void updateProgressStep(int step, String message) {
		currentStep = step;
		double progress = (double) step / totalSteps;
		libraryProgressBar.setProgress(progress);

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
