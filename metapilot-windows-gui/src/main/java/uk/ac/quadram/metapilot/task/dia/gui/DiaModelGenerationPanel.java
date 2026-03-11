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
import javafx.scene.control.cell.CheckBoxTableCell;
import javafx.beans.property.SimpleStringProperty;
import javafx.beans.property.StringProperty;
import javafx.beans.property.SimpleBooleanProperty;
import javafx.beans.property.BooleanProperty;
import javafx.concurrent.Task;
import uk.ac.quadram.metapilot.task.mag.MagDbItem;
import uk.ac.quadram.metapilot.task.mag.MagDbItem.MagModel;
import uk.ac.quadram.metapilot.task.dia.DiaPyTorchModelTask;
import uk.ac.quadram.metapilot.task.common.gui.AbstractWorkflowPanel;
import uk.ac.quadram.metapilot.task.dia.par.MetaParameterDia;
import uk.ac.quadram.metapilot.resources.MetaResourceParameters;
import uk.ac.quadram.metapilot.task.mag.MagDbConfigIO;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.Optional;

import org.apache.logging.log4j.LogManager;

import javafx.animation.KeyFrame;
import javafx.animation.Timeline;
import javafx.util.Duration;

public class DiaModelGenerationPanel extends VBox {

	// Model table components
	private TableView<ModelRow> modelTable;
	private ObservableList<ModelRow> modelData = FXCollections.observableArrayList();

	// Generation components
	private TextField intensityFileField;
	private TextField modelNameField;

	// Current MAG database reference
	private MagDbItem currentMagDbItem;

	// Console output for feedback
	private TextArea consoleOutput;

	// Internal progress tracking
	private ProgressBar modelProgressBar;

	// Callback interface for parameter updates
	public interface ParameterUpdateCallback {
		void onParametersChanged();
	}

	private ParameterUpdateCallback parameterCallback;

	// Required parameters for DiaPyTorchModelTask
	private MetaParameterDia metaParameter;
	private MetaResourceParameters resourceParameters;
	private File configFile;
	private AbstractWorkflowPanel<MetaParameterDia> parentPanel;

	// UI components
	private Button generateModelButton;

	public DiaModelGenerationPanel() {
		initializeUI();
	}

	public DiaModelGenerationPanel(TextArea consoleOutput) {
		this.consoleOutput = consoleOutput;
		initializeUI();
	}

	public DiaModelGenerationPanel(MetaParameterDia metaParameter, MetaResourceParameters resourceParameters,
			File configFile) {
		this.metaParameter = metaParameter;
		this.resourceParameters = resourceParameters;
		this.configFile = configFile;

		// CRITICAL: Initialize MagDbItem from metaParameter immediately
		if (metaParameter != null && metaParameter.getMagDbItem() != null) {
			this.currentMagDbItem = metaParameter.getMagDbItem();
			System.out.println("DEBUG: DiaModelGenerationPanel constructor - initialized MagDbItem: "
					+ currentMagDbItem.getCatalogueID());
		} else {
			System.err.println("ERROR: DiaModelGenerationPanel constructor - MetaParameterDia or MagDbItem is null!");
		}

		initializeUI();

		// Load available models after UI is initialized
		loadAvailableModels();
	}

	public DiaModelGenerationPanel(MetaParameterDia metaParameter, MetaResourceParameters resourceParameters,
			File configFile, TextArea consoleOutput) {
		this.metaParameter = metaParameter;
		this.resourceParameters = resourceParameters;
		this.configFile = configFile;
		this.consoleOutput = consoleOutput;

		// CRITICAL: Initialize MagDbItem from metaParameter immediately
		if (metaParameter != null && metaParameter.getMagDbItem() != null) {
			this.currentMagDbItem = metaParameter.getMagDbItem();
			System.out.println("DEBUG: DiaModelGenerationPanel constructor - initialized MagDbItem: "
					+ currentMagDbItem.getCatalogueID());
		} else {
			System.err.println("ERROR: DiaModelGenerationPanel constructor - MetaParameterDia or MagDbItem is null!");
		}

		initializeUI();

		// Load available models after UI is initialized
		loadAvailableModels();
	}

	public void setParameterUpdateCallback(ParameterUpdateCallback callback) {
		this.parameterCallback = callback;
	}

	public void setConsoleOutput(TextArea consoleOutput) {
		this.consoleOutput = consoleOutput;
	}

	public void setMagDbItem(MagDbItem magDbItem) {
		this.currentMagDbItem = magDbItem;
		loadAvailableModels();
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

	private void initializeUI() {
		setSpacing(10);
		setPadding(new Insets(10));
		setStyle("-fx-background-color: #f0fff0; -fx-background-radius: 5;");
		setMinHeight(400);
		VBox.setVgrow(this, Priority.ALWAYS);

		Label titleLabel = createSectionTitle("ML Model Generation (Optional)");

		// Optional hint
		Label optionalHint = new Label(
				"⚠ Model generation is NOT mandatory - it's an optional feature for enhanced identification");
		optionalHint.setStyle("-fx-font-size: 11px; -fx-text-fill: #ff6600; -fx-font-weight: bold;");
		optionalHint.setWrapText(true);

		// File selection for peptide intensity table
		VBox inputBox = createInputSection();

		// Generation controls
		HBox generationBox = createGenerationControls();

		// Available models table
		VBox tableSection = createModelTable();

		// Progress and console section
		VBox progressSection = createProgressSection();

		getChildren().addAll(titleLabel, optionalHint, inputBox, generationBox, tableSection, progressSection);
	}

	private VBox createInputSection() {
		VBox inputBox = new VBox(8);
		Label inputLabel = new Label("Peptide Intensity Table for Model Generation:");
		inputLabel.setStyle("-fx-font-weight: bold;");

		HBox fileBox = new HBox(10);
		fileBox.setAlignment(Pos.CENTER_LEFT);
		intensityFileField = new TextField();
		intensityFileField.setPromptText("Select peptide intensity table file (.txt, .tsv, .csv)");
		intensityFileField.setPrefWidth(280);
		intensityFileField.setEditable(false);

		Button browseIntensityButton = new Button("Browse");
		browseIntensityButton.setOnAction(e -> browseIntensityFile());

		fileBox.getChildren().addAll(intensityFileField, browseIntensityButton);
		inputBox.getChildren().addAll(inputLabel, fileBox);
		return inputBox;
	}

	private HBox createGenerationControls() {
		HBox generationBox = new HBox(10);
		generationBox.setAlignment(Pos.CENTER_LEFT);

		modelNameField = new TextField();
		modelNameField.setPromptText("Model name");
		modelNameField.setPrefWidth(150);

		generateModelButton = new Button("Generate Model");
		generateModelButton.setPrefWidth(130);
		generateModelButton.setOnAction(e -> generateModel());

		generationBox.getChildren().addAll(new Label("Model Name:"), modelNameField, generateModelButton);
		return generationBox;
	}

	private VBox createModelTable() {
		VBox tableSection = new VBox(10);

		// Header with controls
		HBox headerBox = new HBox(10);
		headerBox.setAlignment(Pos.CENTER_LEFT);
		Label availableLabel = new Label("Available Models (select for identification - optional):");
		availableLabel.setStyle("-fx-font-weight: bold;");

		Button checkAvailableButton = new Button("Check Available");
		checkAvailableButton.setTooltip(new Tooltip("Auto-detect models in MAG database folder"));
		checkAvailableButton.setOnAction(e -> checkAvailableModels());

		Button deleteButton = new Button("Delete Selected");
		deleteButton.setOnAction(e -> deleteSelectedModel());

		headerBox.getChildren().addAll(availableLabel, checkAvailableButton, deleteButton);

		// Table
		modelTable = new TableView<>(modelData);
		modelTable.setMinHeight(120);
		modelTable.setPrefHeight(150);
		modelTable.setMaxHeight(Double.MAX_VALUE);
		modelTable.setEditable(true);
		modelTable.setColumnResizePolicy(TableView.UNCONSTRAINED_RESIZE_POLICY);
		VBox.setVgrow(modelTable, Priority.ALWAYS);

		// Use column (checkbox)
		TableColumn<ModelRow, Boolean> useCol = new TableColumn<>("Use");
		useCol.setCellValueFactory(new PropertyValueFactory<>("use"));
		useCol.setCellFactory(CheckBoxTableCell.forTableColumn(useCol));
		useCol.setOnEditCommit(e -> notifyParameterUpdate());
		useCol.setMinWidth(40);
		useCol.setPrefWidth(50);

		// Model Name column
		TableColumn<ModelRow, String> nameCol = new TableColumn<>("Model Name");
		nameCol.setCellValueFactory(new PropertyValueFactory<>("modelName"));
		nameCol.setMinWidth(120);
		nameCol.setPrefWidth(150);

		// Training Size column
		TableColumn<ModelRow, String> trainSizeCol = new TableColumn<>("Train Size");
		trainSizeCol.setCellValueFactory(new PropertyValueFactory<>("trainSize"));
		trainSizeCol.setMinWidth(80);
		trainSizeCol.setPrefWidth(90);

		// Model Size column
		TableColumn<ModelRow, String> modelSizeCol = new TableColumn<>("Model Size");
		modelSizeCol.setCellValueFactory(new PropertyValueFactory<>("modelSize"));
		modelSizeCol.setMinWidth(80);
		modelSizeCol.setPrefWidth(90);

		// Correlation column
		TableColumn<ModelRow, String> corrCol = new TableColumn<>("Correlation");
		corrCol.setCellValueFactory(new PropertyValueFactory<>("corr"));
		corrCol.setMinWidth(80);
		corrCol.setPrefWidth(90);

		// Median column
		TableColumn<ModelRow, String> medCol = new TableColumn<>("Median");
		medCol.setCellValueFactory(new PropertyValueFactory<>("med"));
		medCol.setMinWidth(80);
		medCol.setPrefWidth(90);

		@SuppressWarnings("unchecked")
		TableColumn<ModelRow, ?>[] columns = new TableColumn[] { useCol, nameCol, trainSizeCol, modelSizeCol, corrCol,
				medCol };
		modelTable.getColumns().addAll(columns);

		tableSection.getChildren().addAll(headerBox, modelTable);
		return tableSection;
	}

	private VBox createProgressSection() {
		VBox progressSection = new VBox(8);
		progressSection.setPadding(new Insets(10));
		progressSection.setStyle("-fx-background-color: #f8f8f8; -fx-background-radius: 5;");
		progressSection.setMinHeight(200); // Increased minimum height
		progressSection.setPrefHeight(400); // Set preferred height to make it bigger
		VBox.setVgrow(progressSection, Priority.ALWAYS); // Allow it to grow and fill space

		modelProgressBar = new ProgressBar(0);
		modelProgressBar.setMaxWidth(Double.MAX_VALUE);
		modelProgressBar.setPrefHeight(25); // Made progress bar taller
		modelProgressBar.setStyle("-fx-accent: #4CAF50;"); // Green color for progress

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

		progressSection.getChildren().addAll(modelProgressBar, consoleTitle, consoleOutput);

		return progressSection;
	}

	// Model Row class based on MagModel structure
	public static class ModelRow {
		private final SimpleBooleanProperty use;
		private final SimpleStringProperty modelName;
		private final SimpleStringProperty trainSize;
		private final SimpleStringProperty modelSize;
		private final SimpleStringProperty corr;
		private final SimpleStringProperty med;
		private MagModel magModel; // Reference to original MagModel

		public ModelRow(boolean use, MagModel magModel) {
			this.use = new SimpleBooleanProperty(use);
			this.modelName = new SimpleStringProperty(magModel.getModelName());
			this.trainSize = new SimpleStringProperty(String.valueOf(magModel.getTrainSize()));
			this.modelSize = new SimpleStringProperty(String.valueOf(magModel.getModelSize()));
			this.corr = new SimpleStringProperty(String.valueOf(magModel.getCorr()));
			this.med = new SimpleStringProperty(String.valueOf(magModel.getMed()));
			this.magModel = magModel;
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

		public String getModelName() {
			return modelName.get();
		}

		public StringProperty modelNameProperty() {
			return modelName;
		}

		public String getTrainSize() {
			return trainSize.get();
		}

		public StringProperty trainSizeProperty() {
			return trainSize;
		}

		public String getModelSize() {
			return modelSize.get();
		}

		public StringProperty modelSizeProperty() {
			return modelSize;
		}

		public String getCorr() {
			return corr.get();
		}

		public StringProperty corrProperty() {
			return corr;
		}

		public String getMed() {
			return med.get();
		}

		public StringProperty medProperty() {
			return med;
		}

		public MagModel getMagModel() {
			return magModel;
		}
	}

	private void browseIntensityFile() {
		javafx.stage.FileChooser chooser = new javafx.stage.FileChooser();
		chooser.setTitle("Select Peptide Intensity Table File");
		chooser.getExtensionFilters()
				.add(new javafx.stage.FileChooser.ExtensionFilter("Intensity Table Files", "*.txt", "*.tsv", "*.csv"));
		File file = chooser.showOpenDialog(getScene().getWindow());
		if (file != null) {
			intensityFileField.setText(file.getAbsolutePath());
		}
	}

	private void generateModel() {
		String intensityFile = intensityFileField.getText();
		String modelName = modelNameField.getText();

		if (intensityFile == null || intensityFile.trim().isEmpty()) {
			showError("Input Required", "Please select a peptide intensity table file.");
			return;
		}

		if (modelName == null || modelName.trim().isEmpty()) {
			showError("Name Required", "Please enter a model name.");
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

		// If currentMagDbItem is still null, try to get from metaParameter
		if (currentMagDbItem == null && metaParameter != null) {
			currentMagDbItem = metaParameter.getMagDbItem();
		}

		// Create and run the real DiaPyTorchModelTask (pass modelName so it's used when
		// registering)
		DiaPyTorchModelTask modelTask = new DiaPyTorchModelTask(metaParameter, resourceParameters, modelProgressBar,
				LogManager.getLogger(DiaModelGenerationPanel.class), intensityFile, modelName, "PyTorch");

		// Set up task event handlers
		modelTask.setOnRunning(e -> {
			Platform.runLater(() -> {
				// Disable the button to prevent multiple clicks
				if (generateModelButton != null) {
					generateModelButton.setDisable(true);
					generateModelButton.setText("Generating...");
				}

				if (consoleOutput != null) {
					consoleOutput.appendText("=== Model Generation Started ===\n");
					consoleOutput.appendText("Model name: " + modelName + "\n");
					consoleOutput.appendText("Input file: " + intensityFile + "\n");
					consoleOutput.appendText("Using PyTorch model for peptide intensity prediction...\n");
					consoleOutput.appendText("Python executable: " + (resourceParameters != null
							? resourceParameters.getToolPath(
									uk.ac.quadram.metapilot.resources.MetaResourceParameters.ToolType.PYTHON)
							: "Not configured") + "\n");
				}
			});
		});

		modelTask.setOnSucceeded(e -> {
			Platform.runLater(() -> {
				// Re-enable the button
				if (generateModelButton != null) {
					generateModelButton.setDisable(false);
					generateModelButton.setText("Generate Model");
				}

				modelProgressBar.setProgress(1.0);

				if (consoleOutput != null) {
					consoleOutput.appendText("Model generation completed successfully: " + modelName + "\n");
					consoleOutput.appendText("Model has been trained and saved to the MAG database.\n");
				}

				// Update the configuration file with the new model information
				if (configFile != null && currentMagDbItem != null) {
					try {
						// Determine if this is a MGnify or local MAG database
						boolean isMgnify = currentMagDbItem.getCatalogueID().startsWith("MGYG")
								|| currentMagDbItem.getCatalogueID().contains("mgnify")
								|| currentMagDbItem.getCatalogueID().contains("MGnify");

						if (consoleOutput != null) {
							consoleOutput
									.appendText("DEBUG: Updating config file: " + configFile.getAbsolutePath() + "\n");
							consoleOutput.appendText("DEBUG: MAG DB ID: " + currentMagDbItem.getCatalogueID() + "\n");
							consoleOutput.appendText("DEBUG: Is MGnify: " + isMgnify + "\n");
							consoleOutput.appendText(
									"DEBUG: Models count: " + currentMagDbItem.getMagModels().length + "\n");
						}

						// Update the configuration file with the updated MagDbItem
						MagDbConfigIO.updateMagDbItem(currentMagDbItem, isMgnify, configFile);

						if (consoleOutput != null) {
							consoleOutput.appendText("✓ Configuration file updated successfully.\n");
						}

					} catch (Exception ex) {
						if (consoleOutput != null) {
							consoleOutput
									.appendText("✗ Failed to update configuration file: " + ex.getMessage() + "\n");
						}
						System.err.println("Failed to update MAG database configuration: " + ex.getMessage());
						ex.printStackTrace();
					}
				} else {
					if (consoleOutput != null) {
						consoleOutput.appendText(
								"✗ Configuration file or MAG database item not set. Configuration not updated.\n");
						consoleOutput.appendText("  Config file: "
								+ (configFile != null ? configFile.getAbsolutePath() : "null") + "\n");
						consoleOutput.appendText("  MAG DB item: "
								+ (currentMagDbItem != null ? currentMagDbItem.getCatalogueID() : "null") + "\n");
					}
				}

				// Refresh the model table
				loadAvailableModels();
				notifyParameterUpdate();

				// Clear input fields
				intensityFileField.clear();
				modelNameField.clear();

				// Reset progress after a delay
				Timeline resetTimeline = new Timeline(new KeyFrame(Duration.seconds(3), ev -> {
					modelProgressBar.setProgress(0);
				}));
				resetTimeline.play();
			});
		});

		modelTask.setOnFailed(e -> {
			Platform.runLater(() -> {
				// Re-enable the button
				if (generateModelButton != null) {
					generateModelButton.setDisable(false);
					generateModelButton.setText("Generate Model");
				}

				modelProgressBar.setProgress(0);

				Throwable exception = modelTask.getException();
				String errorMessage = exception != null ? exception.getMessage() : "Unknown error occurred";

				if (consoleOutput != null) {
					consoleOutput.appendText("=== Model Generation Failed ===\n");
					consoleOutput.appendText("✗ Error: " + errorMessage + "\n");
					consoleOutput.appendText("✗ Check the console output above for Python error details.\n");
					consoleOutput.appendText("✗ Common issues:\n");
					consoleOutput.appendText("  - Missing Python modules (torch, pandas, numpy, etc.)\n");
					consoleOutput.appendText("  - Incorrect Python path configuration\n");
					consoleOutput.appendText("  - Insufficient disk space or permissions\n");
				}

				showError("Model Generation Failed", "Failed to generate model: " + errorMessage
						+ "\n\nCheck the console output for detailed Python error information.");
			});
		});

		modelTask.setOnCancelled(e -> {
			Platform.runLater(() -> {
				// Re-enable the button
				if (generateModelButton != null) {
					generateModelButton.setDisable(false);
					generateModelButton.setText("Generate Model");
				}

				modelProgressBar.setProgress(0);

				if (consoleOutput != null) {
					consoleOutput.appendText("Model generation was cancelled.\n");
				}
			});
		});

		// Set up progress monitoring from the task
		modelTask.messageProperty().addListener((obs, oldMsg, newMsg) -> {
			if (newMsg != null && !newMsg.isEmpty()) {
				Platform.runLater(() -> {
					if (consoleOutput != null) {
						consoleOutput.appendText("Progress: " + newMsg + "\n");
					}
				});
			}
		});

		// Start the task
		new Thread(modelTask).start();
	}

	private void checkAvailableModels() {
		if (currentMagDbItem == null) {
			showError("No Database", "Please select a MAG database first.");
			return;
		}

		Task<Void> checkTask = new Task<Void>() {
			@Override
			protected Void call() throws Exception {
				Platform.runLater(() -> {
					if (consoleOutput != null) {
						consoleOutput.appendText("Scanning for available models in MAG database folder...\n");
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

					// Check the models folder
					File modelsFolder = new File(magDbFolder, "models");
					if (!modelsFolder.exists()) {
						Platform.runLater(() -> {
							if (consoleOutput != null) {
								consoleOutput.appendText(
										"No models folder found at: " + modelsFolder.getAbsolutePath() + "\n");
							}
						});
						return null;
					}

					Platform.runLater(() -> {
						if (consoleOutput != null) {
							consoleOutput
									.appendText("Scanning models folder: " + modelsFolder.getAbsolutePath() + "\n");
						}
					});

					// Scan for model folders
					File[] modelFolders = modelsFolder.listFiles(File::isDirectory);
					if (modelFolders == null || modelFolders.length == 0) {
						Platform.runLater(() -> {
							if (consoleOutput != null) {
								consoleOutput.appendText("No model folders found in models directory.\n");
							}
						});
						return null;
					}

					final int[] foundModels = { 0 };
					for (File modelFolder : modelFolders) {
						String modelName = modelFolder.getName();

						// Look for _pep_rank.tsv file in the model folder (this is the key output file)
						File pepRankFile = new File(modelFolder, modelName + "_pep_rank.tsv");
						File oldPepRankFile = new File(modelFolder, modelName + "pep_rank.tsv"); // Old naming without
																									// underscore
						File namedPytorchPredictionsFile = new File(modelFolder,
								modelName + "_pytorch_predictions.tsv");
						File oldPytorchPredictionsFile = new File(modelFolder, "pytorch_predictions.tsv");
						File trainDataFile = new File(modelFolder, "trainData.tsv");
						File predictionDataFile = new File(modelFolder, "predictionData.tsv");

						// Check if the key model files exist (support both old and new naming)
						boolean hasRankFile = pepRankFile.exists() || oldPepRankFile.exists();
						boolean hasPredictionsFile = namedPytorchPredictionsFile.exists()
								|| oldPytorchPredictionsFile.exists();

						if (hasRankFile && hasPredictionsFile) {
							Platform.runLater(() -> {
								if (consoleOutput != null) {
									consoleOutput.appendText("Found model: " + modelName + " (pep_rank: "
											+ pepRankFile.exists() + ", old_pep_rank: " + oldPepRankFile.exists()
											+ ", named_predictions: " + namedPytorchPredictionsFile.exists()
											+ ", old_predictions: " + oldPytorchPredictionsFile.exists()
											+ ", trainData: " + trainDataFile.exists() + ", predictionData: "
											+ predictionDataFile.exists() + ")\n");
								}
							});

							// Add the model to the MagDbItem
							System.out.println(
									"DEBUG: Attempting to add model from folder: " + modelFolder.getAbsolutePath());

							// Check if the required info file exists
							File infoFile = new File(modelFolder, modelName + "_info.txt");
							System.out.println("DEBUG: Looking for info file: " + infoFile.getAbsolutePath()
									+ " (exists: " + infoFile.exists() + ")");

							if (infoFile.exists()) {
								currentMagDbItem.addModelFromFile(modelFolder);
								System.out.println("DEBUG: Successfully added model via addModelFromFile()");
							} else {
								// Try to read from pytorch_model_results.txt instead
								File resultsFile = new File(modelFolder, "pytorch_model_results.txt");
								if (resultsFile.exists()) {
									try {
										// Read actual values from pytorch_model_results.txt
										int trainSize = 0;
										double corr = 0.0;
										double med = 0.0; // Keep as placeholder

										try (BufferedReader reader = new BufferedReader(new FileReader(resultsFile))) {
											String line;
											boolean inDatasetSizes = false;
											while ((line = reader.readLine()) != null) {
												line = line.trim();

												// Check for Dataset Sizes section
												if (line.equals("Dataset Sizes:")) {
													inDatasetSizes = true;
													continue;
												}

												// Parse Training samples
												if (inDatasetSizes && line.startsWith("Training samples:")) {
													try {
														trainSize = Integer.parseInt(
																line.substring("Training samples:".length()).trim());
													} catch (NumberFormatException e) {
														System.out.println(
																"DEBUG: Could not parse training samples from: "
																		+ line);
													}
												}

												// Parse Correlation from metrics section
												if (line.startsWith("Correlation:")) {
													try {
														corr = Double.parseDouble(
																line.substring("Correlation:".length()).trim());
													} catch (NumberFormatException e) {
														System.out.println(
																"DEBUG: Could not parse correlation from: " + line);
													}
												}
											}
										}

										// Fallback to counting lines if parsing failed
										if (trainSize == 0) {
											trainSize = countLines(new File(modelFolder, "trainData.tsv")) - 1;
										}

										int dataSize = countLines(new File(modelFolder, "predictionData.tsv")) - 1;

										// Add model manually with parsed values
										currentMagDbItem.addModel(modelName, trainSize, dataSize, corr, med);

										System.out.println("DEBUG: Added model from results file - Name: " + modelName
												+ ", TrainSize: " + trainSize + ", DataSize: " + dataSize
												+ ", Correlation: " + corr);

									} catch (Exception e) {
										System.out.println("DEBUG: Error reading results file for model " + modelName
												+ ": " + e.getMessage());
										// Fallback to placeholder values
										currentMagDbItem.addModel(modelName, 0, 0, 0.0, 0.0);
									}
								} else {
									// Create a basic model entry without the info file
									// We'll use placeholder values since we don't have the info file
									System.out.println(
											"DEBUG: Neither info file nor results file found, creating basic model entry");
									currentMagDbItem.addModel(modelName, 0, 0, 0.0, 0.0);
								}
							}

							foundModels[0]++;
						} else {
							Platform.runLater(() -> {
								if (consoleOutput != null) {
									consoleOutput.appendText("Skipping folder (no model files): " + modelName + "\n");
								}
							});
						}
					}

					Platform.runLater(() -> {
						if (consoleOutput != null) {
							consoleOutput.appendText("Model scan completed. Found " + foundModels[0] + " models.\n");
						}
						loadAvailableModels();
					});

				} catch (Exception e) {
					Platform.runLater(() -> {
						if (consoleOutput != null) {
							consoleOutput.appendText("ERROR during model scan: " + e.getMessage() + "\n");
						}
					});
					System.err.println("Error scanning models: " + e.getMessage());
					e.printStackTrace();
				}

				return null;
			}
		};

		new Thread(checkTask).start();
	}

	private void deleteSelectedModel() {
		ModelRow selected = modelTable.getSelectionModel().getSelectedItem();
		if (selected == null) {
			showError("No Selection", "Please select a model to delete.");
			return;
		}

		// Create confirmation dialog with three options
		Alert alert = new Alert(Alert.AlertType.CONFIRMATION);
		alert.setTitle("Delete Model");
		alert.setHeaderText("Delete model: " + selected.getModelName());
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
					currentMagDbItem.removeModel(selected.getMagModel());
				}
				modelData.remove(selected);
				notifyParameterUpdate();

				if (consoleOutput != null) {
					consoleOutput.appendText("Model removed from table: " + selected.getModelName() + "\n");
				}
			} else if (result.get() == deleteFilesButton) {
				// Remove from table, MAG database, and delete files
				if (currentMagDbItem != null) {
					currentMagDbItem.removeModel(selected.getMagModel());
				}
				// Note: MagModel doesn't have a delete() method like MagPepLib,
				// but in real implementation, you would add file deletion logic here
				modelData.remove(selected);
				notifyParameterUpdate();

				if (consoleOutput != null) {
					consoleOutput.appendText("Model and files deleted: " + selected.getModelName() + "\n");
				}
			}
			// Cancel - do nothing
		}
	}

	private void loadAvailableModels() {
		modelData.clear();
		if (currentMagDbItem != null) {
			MagModel[] magModels = currentMagDbItem.getMagModels();
			System.out.println("DEBUG: Loading models from MagDbItem " + currentMagDbItem.getCatalogueID() + " - found "
					+ magModels.length + " models");

			for (MagModel model : magModels) {
				modelData.add(new ModelRow(false, model));
				System.out.println("DEBUG: Added model to table: " + model.getModelName() + " (train size: "
						+ model.getTrainSize() + ")");
			}

			if (consoleOutput != null) {
				consoleOutput.appendText("Loaded " + magModels.length + " models from MAG database.\n");
			}
		} else {
			System.out.println("DEBUG: Cannot load models - currentMagDbItem is null");
			if (consoleOutput != null) {
				consoleOutput.appendText("No MAG database selected - cannot load models.\n");
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

	private int countLines(File f) {
		int count = 0;
		if (f == null || !f.exists())
			return 0;
		try (BufferedReader br = new BufferedReader(new FileReader(f))) {
			while (br.readLine() != null)
				count++;
		} catch (IOException e) {
			return 0;
		}
		return count;
	}

	private void notifyParameterUpdate() {
		System.out.println("DEBUG: DiaModelGenerationPanel.notifyParameterUpdate() called");
		if (parameterCallback != null) {
			System.out.println("DEBUG: Calling parameterCallback.onParametersChanged()");
			parameterCallback.onParametersChanged();
		} else {
			System.out.println("DEBUG: parameterCallback is null");
		}
	}

	// Public methods for accessing selected models
	public MagModel[] getSelectedModels() {
		return modelData.stream().filter(ModelRow::getUse).map(ModelRow::getMagModel).toArray(MagModel[]::new);
	}

	public String[] getSelectedModelNames() {
		return modelData.stream().filter(ModelRow::getUse).map(ModelRow::getModelName).toArray(String[]::new);
	}
}
