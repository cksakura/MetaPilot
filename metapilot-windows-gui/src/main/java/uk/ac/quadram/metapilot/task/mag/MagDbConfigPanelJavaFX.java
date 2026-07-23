package uk.ac.quadram.metapilot.task.mag;

import javafx.application.Platform;
import javafx.beans.property.SimpleBooleanProperty;
import javafx.beans.property.SimpleIntegerProperty;
import javafx.beans.property.SimpleStringProperty;
import javafx.concurrent.Task;
import javafx.event.ActionEvent;
import javafx.geometry.Insets;
import javafx.geometry.Pos;
import javafx.scene.Node;
import javafx.scene.control.*;
import javafx.scene.control.cell.TextFieldTableCell;
import javafx.scene.layout.*;
import javafx.stage.DirectoryChooser;
import javafx.scene.text.Font;
import javafx.scene.text.FontWeight;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Date;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import uk.ac.quadram.metapilot.core.MGYG.HapExtractor;
import uk.ac.quadram.metapilot.core.MGYG.MGYGFtpTaskJavaFX;

/**
 * JavaFX implementation of the MagDbConfigPanel Provides the same functionality
 * as the Swing version but with JavaFX components
 */
public class MagDbConfigPanelJavaFX extends BorderPane {

	private MagDbItem[][] magDbItems;
	private TableView<MagDbItem> mgnifyTable;
	private TableView<MagDbItem> localTable;

	// MGnify panel controls
	private ComboBox<MagDbItem> mgnifyComboBox;
	private TextField spCountTextField;
	private TextField latestVerTextField;
	private TextField mgnifyOldVerTextField;
	private TextField mgnifyRepoTextField;
	private Button downloadButton;
	private ProgressBar mgnifyDLBar;
	private TextArea mgnifyConsoleArea;

	// Local panel controls
	private TextField localNameTextField;
	private TextField identifierTextField;
	private TextField localRepoTextField;
	private TextField localVersionTextField;
	private TextField localDbTextField;
	private TextField localFuncTextField;
	private TextField localTaxaTextField;

	private SimpleDateFormat format = new SimpleDateFormat("MM_dd_yyyy_HH_mm_ss");

	/**
	 * Create the panel with the given database items
	 */
	public MagDbConfigPanelJavaFX(MagDbItem[][] magDbItems) {
		this.magDbItems = magDbItems;
		setPadding(new Insets(10));

		// Create the main content panels
		TabPane tabPane = new TabPane();

		// MGnify Tab
		Tab mgnifyTab = new Tab("MGnify Database");
		mgnifyTab.setClosable(false);
		mgnifyTab.setContent(createMgnifyPanel());

		// Local Tab
		Tab localTab = new Tab("Local MAG Database");
		localTab.setClosable(false);
		localTab.setContent(createLocalPanel());

		tabPane.getTabs().addAll(mgnifyTab, localTab);

		// Set the tab pane as the center content
		setCenter(tabPane);
	}

	/**
	 * Creates the MGnify database panel
	 */
	private Pane createMgnifyPanel() {
		GridPane gridPane = new GridPane();
		gridPane.setHgap(20);
		gridPane.setVgap(10);
		gridPane.setPadding(new Insets(15));

		// Left side: MGnify database list
		VBox leftPanel = createMGnifyListPanel();

		// Right side: MGnify download panel
		VBox rightPanel = createMGnifyDownloadPanel();

		// Add panels to grid
		gridPane.add(leftPanel, 0, 0);
		gridPane.add(rightPanel, 1, 0);

		// Allow panels to grow to fill available space
		GridPane.setHgrow(leftPanel, Priority.ALWAYS);
		GridPane.setHgrow(rightPanel, Priority.ALWAYS);
		GridPane.setVgrow(leftPanel, Priority.ALWAYS);
		GridPane.setVgrow(rightPanel, Priority.ALWAYS);

		return gridPane;
	}

	/**
	 * Creates the MGnify database list panel
	 */
	private VBox createMGnifyListPanel() {
		VBox panel = new VBox(10);
		panel.setPadding(new Insets(10));
		panel.setStyle("-fx-border-color: #cccccc; -fx-border-radius: 5;");

		Label titleLabel = new Label("MGnify Database List");
		titleLabel.setFont(Font.font("System", FontWeight.BOLD, 14));

		// Create table view for MGnify databases
		mgnifyTable = new TableView<>();
		mgnifyTable.setEditable(false);
		mgnifyTable.setPlaceholder(new Label("No MGnify databases available"));

		// Create columns
		TableColumn<MagDbItem, Boolean> availableColumn = new TableColumn<>("Available");
		availableColumn.setCellValueFactory(data -> new SimpleBooleanProperty(data.getValue().isAvailable()));
		availableColumn.setCellFactory(col -> new TableCell<MagDbItem, Boolean>() {
			private final CheckBox checkBox = new CheckBox();

			{
				// Make checkbox non-editable
				checkBox.setDisable(true);
				// Center the checkbox in the cell
				setAlignment(Pos.CENTER);
			}

			@Override
			protected void updateItem(Boolean item, boolean empty) {
				super.updateItem(item, empty);
				if (empty || item == null) {
					setGraphic(null);
				} else {
					checkBox.setSelected(item);
					setGraphic(checkBox);
				}
			}
		});
		availableColumn.setPrefWidth(80);

		TableColumn<MagDbItem, String> idColumn = new TableColumn<>("Catalogue ID");
		idColumn.setCellValueFactory(data -> new SimpleStringProperty(data.getValue().getCatalogueID()));
		idColumn.setPrefWidth(120);

		TableColumn<MagDbItem, String> nameColumn = new TableColumn<>("Identifier");
		nameColumn.setCellValueFactory(data -> new SimpleStringProperty(data.getValue().getIdentifier()));
		nameColumn.setPrefWidth(150);

		TableColumn<MagDbItem, Number> speciesColumn = new TableColumn<>("Species Count");
		speciesColumn.setCellValueFactory(data -> new SimpleIntegerProperty(data.getValue().getSpeciesCount()));
		speciesColumn.setPrefWidth(100);

		// Add columns to table
		mgnifyTable.getColumns().addAll(availableColumn, idColumn, nameColumn, speciesColumn);

		// Add data to table
		loadMgnifyTableData();

		// Set selection handler
		mgnifyTable.getSelectionModel().selectedItemProperty().addListener((obs, oldVal, newVal) -> {
			if (newVal != null) {
				updateMgnifySelection(newVal);
			}
		});

		// Add Check MGnify button and progress bar
		HBox buttonBox = new HBox(10);
		Button checkButton = new Button("Check MGnify");
		ProgressBar progressBar = new ProgressBar(0);
		progressBar.setPrefWidth(300);
		buttonBox.getChildren().addAll(checkButton, progressBar);

		// Set up the Check MGnify action
		checkButton.setOnAction(e -> {
			checkButton.setDisable(true);
			progressBar.setProgress(ProgressBar.INDETERMINATE_PROGRESS);

			// Create MGYGFtpTaskJavaFX for repository discovery
			MGYGFtpTaskJavaFX checkTask = new MGYGFtpTaskJavaFX();

			// Bind progress bar to task progress
			progressBar.progressProperty().bind(checkTask.progressProperty());

			// Handle task completion
			checkTask.setOnSucceeded(event -> {
				checkButton.setDisable(false);
				progressBar.progressProperty().unbind();
				progressBar.setProgress(0);

				// Update the table with discovered repositories
				MagDbItem[] repositories = checkTask.getRepositoryItems();
				if (repositories != null && repositories.length > 0) {
					// Update magDbItems array
					if (magDbItems == null) {
						magDbItems = new MagDbItem[2][];
					}
					magDbItems[0] = repositories;

					// Refresh the table
					Platform.runLater(() -> {
						loadMgnifyTableData();
						loadMgnifyComboBox();
					});

					showAlert(Alert.AlertType.INFORMATION, "Check Complete",
							"Found " + repositories.length + " MGnify databases");
				} else {
					showAlert(Alert.AlertType.WARNING, "Check Complete", "No MGnify databases found");
				}
			});

			checkTask.setOnFailed(event -> {
				checkButton.setDisable(false);
				progressBar.progressProperty().unbind();
				progressBar.setProgress(0);

				Throwable exception = checkTask.getException();
				showAlert(Alert.AlertType.ERROR, "Check Failed", "Failed to check MGnify databases: "
						+ (exception != null ? exception.getMessage() : "Unknown error"));
			});

			checkTask.setOnCancelled(event -> {
				checkButton.setDisable(false);
				progressBar.progressProperty().unbind();
				progressBar.setProgress(0);
				showAlert(Alert.AlertType.INFORMATION, "Check Cancelled", "MGnify database check was cancelled");
			});

			// Start the task
			Thread thread = new Thread(checkTask);
			thread.setDaemon(true);
			thread.start();
		});

		// Make table take up all available space
		VBox.setVgrow(mgnifyTable, Priority.ALWAYS);

		// Add all components to panel
		panel.getChildren().addAll(titleLabel, mgnifyTable, buttonBox);

		return panel;
	}

	/**
	 * Loads data into the MGnify table
	 */
	private void loadMgnifyTableData() {
		mgnifyTable.getItems().clear();

		if (magDbItems != null && magDbItems[0] != null) {
			for (MagDbItem item : magDbItems[0]) {
				mgnifyTable.getItems().add(item);
			}
		}

		// If any items in the table, select the first one
		if (!mgnifyTable.getItems().isEmpty()) {
			mgnifyTable.getSelectionModel().select(0);
			updateMgnifySelection(mgnifyTable.getItems().get(0));
		}
	}

	/**
	 * Creates the MGnify download panel
	 */
	private VBox createMGnifyDownloadPanel() {
		VBox panel = new VBox(10);
		panel.setPadding(new Insets(10));
		panel.setStyle("-fx-border-color: #cccccc; -fx-border-radius: 5;");

		Label titleLabel = new Label("Download MGnify Database");
		titleLabel.setFont(Font.font("System", FontWeight.BOLD, 14));

		// Create the form grid
		GridPane formGrid = new GridPane();
		formGrid.setHgap(10);
		formGrid.setVgap(10);
		formGrid.setPadding(new Insets(5));

		// Database selection
		Label nameLabel = new Label("Name:");
		mgnifyComboBox = new ComboBox<>();
		mgnifyComboBox.setPrefWidth(200);
		loadMgnifyComboBox();

		// Species count
		Label speciesLabel = new Label("Species count:");
		spCountTextField = new TextField();
		spCountTextField.setEditable(false);

		// Latest version
		Label latestVerLabel = new Label("Latest version:");
		latestVerTextField = new TextField();
		latestVerTextField.setEditable(false);

		// Local version
		Label localVerLabel = new Label("Local version:");
		mgnifyOldVerTextField = new TextField();
		mgnifyOldVerTextField.setEditable(false);

		// Repository path
		Label repoLabel = new Label("Local repository:");
		mgnifyRepoTextField = new TextField();
		Button browseButton = new Button("Browse");

		// Set up browse button action
		browseButton.setOnAction(e -> {
			DirectoryChooser dirChooser = new DirectoryChooser();
			dirChooser.setTitle("Select MGnify Repository Directory");

			// If there's already a path, start from there
			if (!mgnifyRepoTextField.getText().isEmpty()) {
				File dir = new File(mgnifyRepoTextField.getText());
				if (dir.getParentFile() != null && dir.getParentFile().exists()) {
					dirChooser.setInitialDirectory(dir.getParentFile());
				}
			}

			File selectedDir = dirChooser.showDialog(getScene().getWindow());
			if (selectedDir != null) {
				MagDbItem item = mgnifyComboBox.getValue();
				if (item != null) {
					// Check if directory name matches catalogue ID
					String fileName = selectedDir.getName();
					if (!fileName.equals(item.getCatalogueID())) {
						selectedDir = new File(selectedDir, item.getCatalogueID());
					}

					mgnifyRepoTextField.setText(selectedDir.getAbsolutePath());
					item.setRepositoryFile(selectedDir);

					// Update download button enablement
					updateDownloadButton();
				}
			}
		});

		// Download button
		downloadButton = new Button("Download/Update");
		downloadButton.setDisable(true);
		mgnifyDLBar = new ProgressBar(0);
		mgnifyDLBar.setPrefWidth(350);

		// Console output area
		mgnifyConsoleArea = new TextArea();
		mgnifyConsoleArea.setEditable(false);
		mgnifyConsoleArea.setPrefHeight(200);

		// Set up download button action
		downloadButton.setOnAction(e -> {
			MagDbItem selectedItem = mgnifyComboBox.getValue();
			if (selectedItem == null) {
				showAlert(Alert.AlertType.WARNING, "Selection Required", "Please select a MGnify database to download");
				return;
			}

			File repoFile = selectedItem.getRepositoryFile();
			if (repoFile == null || !repoFile.getParentFile().exists()) {
				showAlert(Alert.AlertType.WARNING, "Repository Required",
						"Please specify a valid local repository path");
				return;
			}

			// Disable UI during download
			downloadButton.setDisable(true);
			browseButton.setDisable(true);
			mgnifyComboBox.setDisable(true);

			// Create download task using MGYGFtpTaskJavaFX
			MGYGFtpTaskJavaFX downloadTask = new MGYGFtpTaskJavaFX(selectedItem);

			// Bind progress and messages
			mgnifyDLBar.progressProperty().bind(downloadTask.progressProperty());
			mgnifyConsoleArea.textProperty().bind(downloadTask.messageProperty());

			// Handle task completion
			downloadTask.setOnSucceeded(event -> {
				// Re-enable UI
				downloadButton.setDisable(false);
				browseButton.setDisable(false);
				mgnifyComboBox.setDisable(false);

				// Unbind properties
				mgnifyDLBar.progressProperty().unbind();
				mgnifyConsoleArea.textProperty().unbind();

				// Set database as available
				selectedItem.setAvailable();
				selectedItem.setCurrentFile(repoFile);

				// Update table to show database is now available
				Platform.runLater(() -> {
					mgnifyTable.refresh();
				});

				// Show success message
				showAlert(Alert.AlertType.INFORMATION, "Download Complete",
						"The MGnify database was downloaded successfully");
			});

			downloadTask.setOnFailed(event -> {
				// Re-enable UI
				downloadButton.setDisable(false);
				browseButton.setDisable(false);
				mgnifyComboBox.setDisable(false);

				// Unbind properties
				mgnifyDLBar.progressProperty().unbind();
				mgnifyConsoleArea.textProperty().unbind();

				// Show error message
				Throwable exception = downloadTask.getException();
				showAlert(Alert.AlertType.ERROR, "Download Failed", "Failed to download the MGnify database: "
						+ (exception != null ? exception.getMessage() : "Unknown error"));
			});

			downloadTask.setOnCancelled(event -> {
				// Re-enable UI
				downloadButton.setDisable(false);
				browseButton.setDisable(false);
				mgnifyComboBox.setDisable(false);

				// Unbind properties
				mgnifyDLBar.progressProperty().unbind();
				mgnifyConsoleArea.textProperty().unbind();

				showAlert(Alert.AlertType.INFORMATION, "Download Cancelled", "MGnify database download was cancelled");
			});

			// Start the task
			Thread thread = new Thread(downloadTask);
			thread.setDaemon(true);
			thread.start();
		});

		// Set up combobox selection handler
		mgnifyComboBox.setOnAction(e -> {
			MagDbItem item = mgnifyComboBox.getValue();
			if (item != null) {
				updateMgnifySelection(item);
			}
		});

		// Add form controls to grid
		int row = 0;
		formGrid.add(nameLabel, 0, row);
		formGrid.add(mgnifyComboBox, 1, row);
		formGrid.add(speciesLabel, 2, row);
		formGrid.add(spCountTextField, 3, row);

		row++;
		formGrid.add(latestVerLabel, 0, row);
		formGrid.add(latestVerTextField, 1, row);
		formGrid.add(localVerLabel, 2, row);
		formGrid.add(mgnifyOldVerTextField, 3, row);

		row++;
		formGrid.add(repoLabel, 0, row);
		formGrid.add(mgnifyRepoTextField, 1, row, 2, 1);
		formGrid.add(browseButton, 3, row);

		row++;
		formGrid.add(downloadButton, 0, row);
		formGrid.add(mgnifyDLBar, 1, row, 3, 1);

		// Make the console area expand to fill available space
		VBox.setVgrow(mgnifyConsoleArea, Priority.ALWAYS);

		// Add components to panel
		panel.getChildren().addAll(titleLabel, formGrid, mgnifyConsoleArea);

		return panel;
	}

	/**
	 * Update UI when a MGnify database is selected
	 */
	private void updateMgnifySelection(MagDbItem selectedItem) {
		if (selectedItem != null) {
			// Update combobox selection if available
			if (mgnifyComboBox != null && mgnifyComboBox.getValue() != selectedItem) {
				mgnifyComboBox.setValue(selectedItem);
			}

			// Update text fields
			if (spCountTextField != null) {
				spCountTextField.setText(String.valueOf(selectedItem.getSpeciesCount()));
			}

			if (latestVerTextField != null) {
				latestVerTextField.setText(selectedItem.getCatalogueVersion());
			}

			if (mgnifyOldVerTextField != null) {
				mgnifyOldVerTextField.setText(selectedItem.getOldVersion());
			}

			// Set repository path if available
			if (mgnifyRepoTextField != null) {
				File repoFile = selectedItem.getRepositoryFile();
				if (repoFile != null && repoFile.exists()) {
					mgnifyRepoTextField.setText(repoFile.getAbsolutePath());
				} else {
					mgnifyRepoTextField.clear();
				}
			}

			// Update download button state
			updateDownloadButton();
		}
	}

	/**
	 * Update the state of the download button based on selection
	 */
	private void updateDownloadButton() {
		if (downloadButton != null && mgnifyComboBox != null) {
			MagDbItem item = mgnifyComboBox.getValue();
			if (item != null) {
				downloadButton.setDisable(false);
			} else {
				downloadButton.setDisable(true);
			}
		}
	}

	/**
	 * Load the MGnify combobox with database items
	 */
	private void loadMgnifyComboBox() {
		mgnifyComboBox.getItems().clear();

		if (magDbItems != null && magDbItems[0] != null) {
			for (MagDbItem item : magDbItems[0]) {
				mgnifyComboBox.getItems().add(item);
			}

			if (!mgnifyComboBox.getItems().isEmpty()) {
				mgnifyComboBox.setValue(mgnifyComboBox.getItems().get(0));
				updateMgnifySelection(mgnifyComboBox.getValue());
			}
		}
	}

	/**
	 * Creates the Local MAG database panel
	 */
	private Pane createLocalPanel() {
		BorderPane panel = new BorderPane();
		panel.setPadding(new Insets(15));

		// Table view for local databases
		localTable = new TableView<>();
		localTable.setEditable(true);
		localTable.setPlaceholder(new Label("No local MAG databases configured"));

		// Create columns
		TableColumn<MagDbItem, String> idColumn = new TableColumn<>("Catalogue ID");
		idColumn.setCellValueFactory(data -> new SimpleStringProperty(data.getValue().getCatalogueID()));
		idColumn.setCellFactory(TextFieldTableCell.forTableColumn());
		idColumn.setOnEditCommit(e -> {
			MagDbItem item = e.getRowValue();
			item.setCatalogueID(e.getNewValue());
		});
		idColumn.setPrefWidth(120);

		TableColumn<MagDbItem, String> nameColumn = new TableColumn<>("Identifier");
		nameColumn.setCellValueFactory(data -> new SimpleStringProperty(data.getValue().getIdentifier()));
		nameColumn.setCellFactory(TextFieldTableCell.forTableColumn());
		nameColumn.setOnEditCommit(e -> {
			MagDbItem item = e.getRowValue();
			item.setIdentifier(e.getNewValue());
		});
		nameColumn.setPrefWidth(150);

		TableColumn<MagDbItem, Number> genomeCountColumn = new TableColumn<>("Genome Count");
		genomeCountColumn.setCellValueFactory(data -> new SimpleIntegerProperty(data.getValue().getSpeciesCount()));
		genomeCountColumn.setPrefWidth(100);

		TableColumn<MagDbItem, String> pathColumn = new TableColumn<>("Path");
		pathColumn.setCellValueFactory(
				data -> new SimpleStringProperty(data.getValue().getRepositoryFile().getAbsolutePath()));
		pathColumn.setPrefWidth(300);

		// Actions column
		TableColumn<MagDbItem, Void> actionColumn = new TableColumn<>("Actions");
		actionColumn.setCellFactory(col -> new TableCell<MagDbItem, Void>() {
			private final Button deleteButton = new Button("Delete");

			{
				// Set up delete button
				deleteButton.setOnAction(e -> {
					MagDbItem item = getTableView().getItems().get(getIndex());

					Alert alert = new Alert(Alert.AlertType.CONFIRMATION);
					alert.setTitle("Confirm Deletion");
					alert.setHeaderText("Delete Database Record");
					alert.setContentText("Are you sure you want to remove this database record?\n\n" + "Catalogue ID: "
							+ item.getCatalogueID() + "\n" + "Identifier: " + item.getIdentifier() + "\n\n"
							+ "NOTE: This action will only remove the database record from the configuration file.\n"
							+ "The actual database files in the database folder will NOT be affected or deleted.");

					alert.showAndWait().ifPresent(result -> {
						if (result == ButtonType.OK) {
							// Remove from table
							getTableView().getItems().remove(item);

							// Update the magDbItems array
							if (magDbItems != null && magDbItems[1] != null) {
								List<MagDbItem> localList = new ArrayList<>();
								for (MagDbItem dbItem : magDbItems[1]) {
									if (!dbItem.equals(item)) {
										localList.add(dbItem);
									}
								}
								magDbItems[1] = localList.isEmpty() ? new MagDbItem[0]
										: localList.toArray(new MagDbItem[0]);

								// Save to file
								saveMagDbConfig();

								// Refresh table to update genome count column
								getTableView().refresh();
							}
						}
					});
				});
			}

			@Override
			protected void updateItem(Void item, boolean empty) {
				super.updateItem(item, empty);
				setGraphic(empty ? null : deleteButton);
			}
		});
		actionColumn.setPrefWidth(100);

		// Add columns to table
		localTable.getColumns().addAll(idColumn, nameColumn, genomeCountColumn, pathColumn, actionColumn);

		// Load data
		loadLocalTableData();

		// Add buttons at the bottom
		HBox buttonBox = new HBox(10);
		buttonBox.setPadding(new Insets(10, 0, 0, 0));

		Button addButton = new Button("Add New Database");
		addButton.setOnAction(e -> {
			showAddLocalDatabaseDialog();
		});

		buttonBox.getChildren().add(addButton);

		// Set table as center, buttons at bottom
		panel.setCenter(localTable);
		panel.setBottom(buttonBox);

		return panel;
	}

	/**
	 * Load data into the local MAG database table
	 */
	private void loadLocalTableData() {
		localTable.getItems().clear();

		if (magDbItems != null && magDbItems[1] != null) {
			for (MagDbItem item : magDbItems[1]) {
				localTable.getItems().add(item);
			}
		}
	}

	private void showAddLocalDatabaseDialog() {
		Dialog<MagDbItem> dialog = new Dialog<>();
		dialog.setTitle("Add Local MAG Database");
		dialog.setHeaderText("Configure and validate the local MAG database before adding it.");

		ButtonType addButtonType = new ButtonType("Add Database", ButtonBar.ButtonData.OK_DONE);
		dialog.getDialogPane().getButtonTypes().addAll(addButtonType, ButtonType.CANCEL);

		GridPane grid = new GridPane();
		grid.setHgap(10);
		grid.setVgap(10);
		grid.setPadding(new Insets(10));

		if (localNameTextField == null) {
			localNameTextField = new TextField();
		}
		if (identifierTextField == null) {
			identifierTextField = new TextField();
		}
		if (localRepoTextField == null) {
			localRepoTextField = new TextField();
		}
		if (localVersionTextField == null) {
			localVersionTextField = new TextField();
		}
		if (localDbTextField == null) {
			localDbTextField = new TextField("original_db");
		}
		if (localFuncTextField == null) {
			localFuncTextField = new TextField("eggNOG");
		}
		if (localTaxaTextField == null) {
			localTaxaTextField = new TextField("genomes-all_metadata.tsv");
		}

		localNameTextField.setPromptText("e.g. human-gut");
		identifierTextField.setPromptText("e.g. MGYG");
		localRepoTextField.setPromptText("Select the MAG database folder");
		localVersionTextField.setPromptText("e.g. v1.0");
		localDbTextField.setPromptText("DB folder name");
		localFuncTextField.setPromptText("Function annotation folder");
		localTaxaTextField.setPromptText("Taxonomy file name");

		Button browseButton = new Button("Browse");
		browseButton.setOnAction(e -> {
			DirectoryChooser dirChooser = new DirectoryChooser();
			dirChooser.setTitle("Select MAG Database Folder");
			if (!localRepoTextField.getText().trim().isEmpty()) {
				File existingDir = new File(localRepoTextField.getText().trim());
				if (existingDir.exists() && existingDir.isDirectory()) {
					dirChooser.setInitialDirectory(existingDir);
				} else if (existingDir.getParentFile() != null && existingDir.getParentFile().exists()) {
					dirChooser.setInitialDirectory(existingDir.getParentFile());
				}
			}
			File selectedDir = dirChooser.showDialog(getScene().getWindow());
			if (selectedDir != null) {
				localRepoTextField.setText(selectedDir.getAbsolutePath());
			}
		});

		int row = 0;
		grid.add(new Label("Catalogue ID:"), 0, row);
		grid.add(localNameTextField, 1, row, 2, 1);
		row++;
		grid.add(new Label("Identifier:"), 0, row);
		grid.add(identifierTextField, 1, row, 2, 1);
		row++;
		grid.add(new Label("Database folder:"), 0, row);
		grid.add(localRepoTextField, 1, row);
		grid.add(browseButton, 2, row);
		row++;
		grid.add(new Label("Version name:"), 0, row);
		grid.add(localVersionTextField, 1, row, 2, 1);
		row++;
		grid.add(new Label("DB name:"), 0, row);
		grid.add(localDbTextField, 1, row, 2, 1);
		row++;
		grid.add(new Label("Function folder:"), 0, row);
		grid.add(localFuncTextField, 1, row, 2, 1);
		row++;
		grid.add(new Label("Taxonomy file:"), 0, row);
		grid.add(localTaxaTextField, 1, row, 2, 1);

		// Create instructions panel
		VBox instructionsPanel = new VBox(5);
		instructionsPanel.setPadding(new Insets(10));
		instructionsPanel.setStyle("-fx-background-color: #f0f0f0; -fx-border-color: #cccccc; -fx-border-radius: 5;");

		Label instructionsTitle = new Label("Database Requirements:");
		instructionsTitle.setFont(Font.font(null, FontWeight.BOLD, 12));

		TextArea instructionsText = new TextArea();
		instructionsText.setEditable(false);
		instructionsText.setWrapText(true);
		instructionsText.setPrefRowCount(12);
		instructionsText.setPrefHeight(250);
		instructionsText
				.setText("1. File Count: The original_db and eggNOG folders must have the same number of files.\n"
						+ "2. Filename Format:\n"
						+ "   - original_db: Files must be named as 'genomeID.faa' or 'genomeID.fasta'\n"
						+ "   - eggNOG: Files must start with the genomeID\n"
						+ "3. File Pairing: Each genomeID must have a matching pair (one in original_db, one in eggNOG)\n"
						+ "4. FASTA Format: All protein headers must follow the format:\n"
						+ "   >genomeID_proteinID protein/gene name\n"
						+ "   - The genomeID in the header must match the filename\n"
						+ "   - Protein/gene name is mandatory\n"
						+ "5. eggNOG Format: Files must be in eggnogMapper exported format (tab-separated)\n"
						+ "6. Taxonomic information Format: A file contains the taxonomic information named\n "
						+ "   'genomes-all_metadata.tsv' (tab-separated) must be included,\n"
						+ "   at least two columns 'Genome' and 'Lineage' are required in this file, for example\n"
						+ "   ---------------------------------\n"
						+ "   Genome    Lineage\n"
						+ "   MGYG000000001 d__Bacteria;p__Firmicutes_A;c__Clostridia;o__Peptostreptococcales;f__Peptostreptococcaceae;g__GCA-900066495;s__GCA-900066495 sp902362365\n"
						+ "   ---------------------------------\n"
						+ "7. HAP Database: After validation, a HAP (High Abundant Proteins) database will be\n"
						+ "   automatically generated from the original_db folder. This file (Marker4.fasta)\n"
						+ "   must be created successfully for the database to be added.");
		instructionsText.setStyle("-fx-background-color: white;");

		instructionsPanel.getChildren().addAll(instructionsTitle, instructionsText);

		// Create validation panel with progress bar and console
		VBox validationPanel = new VBox(10);
		validationPanel.setPadding(new Insets(10));

		Label validationLabel = new Label("Validation:");
		validationLabel.setPrefWidth(100);
		validationLabel.setFont(Font.font(null, FontWeight.BOLD, 12));

		ProgressBar validationProgressBar = new ProgressBar(0);
		validationProgressBar.setPrefWidth(800);
		validationProgressBar.setMaxWidth(Double.MAX_VALUE);
		validationProgressBar.setMinHeight(18);

		TextArea validationConsole = new TextArea();
		validationConsole.setEditable(false);
		validationConsole.setPrefHeight(200);
		validationConsole.setWrapText(true);

		// Make console scrollable
		ScrollPane consoleScrollPane = new ScrollPane(validationConsole);
		consoleScrollPane.setFitToWidth(true);
		consoleScrollPane.setFitToHeight(true);
		consoleScrollPane.setVbarPolicy(ScrollPane.ScrollBarPolicy.AS_NEEDED);
		consoleScrollPane.setHbarPolicy(ScrollPane.ScrollBarPolicy.NEVER);

		// Create top row with two columns
		HBox topRow = new HBox(10);
		topRow.setPadding(new Insets(10));

		// Left column: Input form (parameters)
		VBox leftColumn = new VBox(10);
		leftColumn.setPrefWidth(450);
		leftColumn.setFillWidth(true);
		leftColumn.getChildren().add(grid);
		HBox.setHgrow(leftColumn, Priority.ALWAYS);
		GridPane.setHgrow(grid, Priority.ALWAYS);

		// Right column: Instructions (requirements)
		VBox rightColumn = new VBox(10);
		rightColumn.setPrefWidth(450);
		rightColumn.getChildren().add(instructionsPanel);
		HBox.setHgrow(rightColumn, Priority.ALWAYS);

		topRow.getChildren().addAll(leftColumn, rightColumn);

		// Bottom section: Progress bar and console
		VBox bottomSection = new VBox(10);
		bottomSection.setPadding(new Insets(10));
		HBox progressRow = new HBox(10, validationLabel, validationProgressBar);
		progressRow.setAlignment(Pos.CENTER_LEFT);
		HBox.setHgrow(validationProgressBar, Priority.ALWAYS);
		bottomSection.getChildren().addAll(progressRow, consoleScrollPane);
		VBox.setVgrow(consoleScrollPane, Priority.ALWAYS);

		// Main layout: Top row and bottom section
		VBox mainContent = new VBox(10);
		mainContent.setPadding(new Insets(10));
		mainContent.getChildren().addAll(topRow, bottomSection);
		VBox.setVgrow(bottomSection, Priority.ALWAYS);

		dialog.getDialogPane().setContent(mainContent);
		dialog.getDialogPane().setPrefSize(900, 700);

		// Store validation result
		final boolean[] validationPassed = { false };
		final MagDbItem[] validatedItem = { null };
		final ValidationTask[] currentValidationTask = { null };

		Node addButton = dialog.getDialogPane().lookupButton(addButtonType);
		addButton.addEventFilter(ActionEvent.ACTION, event -> {
			// Always consume the event first
			event.consume();

			// If validation already passed, allow dialog to close
			if (validationPassed[0] && validatedItem[0] != null) {
				dialog.setResult(validatedItem[0]);
				return;
			}

			// Disable the button during validation
			addButton.setDisable(true);
			validationProgressBar.setProgress(0);
			validationConsole.clear();

			// Run validation in background task
			ValidationTask validationTask = new ValidationTask(localNameTextField.getText().trim(),
					identifierTextField.getText().trim(), localRepoTextField.getText().trim(),
					localDbTextField.getText().trim(), localFuncTextField.getText().trim(),
					localTaxaTextField.getText().trim(), validationConsole, validationProgressBar);

			currentValidationTask[0] = validationTask;

			validationTask.setOnSucceeded(e -> {
				if (validationTask.getValue()) {
					// Validation passed
					validationPassed[0] = true;
					validatedItem[0] = buildLocalMagDbItem(validationTask.getGenomeCount());
					Platform.runLater(() -> {
						// Show success message
						Alert successAlert = new Alert(Alert.AlertType.INFORMATION);
						successAlert.setTitle("Validation Successful");
						successAlert.setHeaderText("Database Validation Passed");
						successAlert.setContentText("All validations passed successfully!\n\n"
								+ "The database has been validated and the HAP database has been generated.\n"
								+ "The database will now be added to the configuration.");
						successAlert.showAndWait();

						// Add the new database immediately and close the dialog
						if (validatedItem[0] != null) {
							applyLocalMagDbItem(validatedItem[0]);
							dialog.close();
						} else {
							addButton.setDisable(false);
						}
					});
				} else {
					// Validation failed, keep dialog open
					addButton.setDisable(false);
					Platform.runLater(() -> {
						Alert failureAlert = new Alert(Alert.AlertType.ERROR);
						failureAlert.setTitle("Validation Failed");
						failureAlert.setHeaderText("Database Validation Failed");
						failureAlert.setContentText("The database validation did not pass.\n\n"
								+ "Please review the validation console output above for detailed error messages.\n\n"
								+ "IMPORTANT: Please read the Database Requirements section at the top of this dialog "
								+ "to understand the required format and structure of the database files.\n\n"
								+ "The database will not be added until all validation checks pass.");
						failureAlert.showAndWait();
					});
				}
			});

			validationTask.setOnFailed(e -> {
				validationConsole.appendText("\nERROR: Validation failed with exception: "
						+ validationTask.getException().getMessage() + "\n");
				addButton.setDisable(false);
				Platform.runLater(() -> {
					Alert failureAlert = new Alert(Alert.AlertType.ERROR);
					failureAlert.setTitle("Validation Error");
					failureAlert.setHeaderText("Validation Process Failed");
					failureAlert.setContentText("An error occurred during validation:\n\n"
							+ validationTask.getException().getMessage() + "\n\n"
							+ "Please review the validation console output above for details.\n\n"
							+ "IMPORTANT: Please read the Database Requirements section at the top of this dialog\n"
							+ "to understand the required format and structure of the database files.");
					failureAlert.showAndWait();
				});
			});

			new Thread(validationTask).start();
		});

		dialog.setResultConverter(dialogButton -> {
			// If validation passed and we have a validated item, return it
			if (dialogButton == addButtonType && validationPassed[0] && validatedItem[0] != null) {
				return validatedItem[0];
			}
			// Otherwise return null (dialog will close with null result)
			return null;
		});

		dialog.showAndWait();
	}

	/**
	 * Validation task that performs comprehensive validation of local MAG database
	 */
	private class ValidationTask extends Task<Boolean> {
		private final String catalogueId;
		private final String identifier;
		private final String dbFolderPath;
		private final String dbName;
		private final String funcName;
		private final String taxaName;
		private final TextArea console;
		private final ProgressBar progressBar;
		private int genomeCount = 0;

		public ValidationTask(String catalogueId, String identifier, String dbFolderPath, String dbName,
				String funcName, String taxaName, TextArea console, ProgressBar progressBar) {
			this.catalogueId = catalogueId;
			this.identifier = identifier;
			this.dbFolderPath = dbFolderPath;
			this.dbName = dbName;
			this.funcName = funcName;
			this.taxaName = taxaName;
			this.console = console;
			this.progressBar = progressBar;
		}

		public int getGenomeCount() {
			return genomeCount;
		}

		@Override
		protected Boolean call() throws Exception {
			Platform.runLater(() -> {
				console.appendText("Starting validation...\n");
				progressBar.setProgress(0.0);
			});

			// Step 1: Basic field validation (10%)
			if (!validateBasicFields()) {
				return false;
			}
			updateProgressBar(0.1, 1.0);

			// Step 1.5: Check for duplicate catalogueId + version (12%)
			if (!validateNoDuplicate()) {
				return false;
			}
			updateProgressBar(0.12, 1.0);

			// Step 2: Folder structure validation (20%)
			File[] folders = validateFolderStructure();
			if (folders == null) {
				return false;
			}
			updateProgressBar(0.2, 1.0);

			// Step 3: File count validation (30%)
			if (!validateFileCounts(folders[0], folders[1])) {
				return false;
			}
			updateProgressBar(0.3, 1.0);

			// Step 4: Filename format and pairing validation (50%)
			Map<String, File[]> pairedFiles = validateFilePairing(folders[0], folders[1]);
			if (pairedFiles == null) {
				return false;
			}
			updateProgressBar(0.5, 1.0);

			// Step 5: FASTA content validation (80%)
			if (!validateFastaContent(pairedFiles)) {
				return false;
			}
			updateProgressBar(0.8, 1.0);

			// Step 6: eggNOG format validation (85%)
			if (!validateEggNogFormat(pairedFiles, folders[1])) {
				return false;
			}
			updateProgressBar(0.85, 1.0);

			// Step 7: Generate HAP database (95%)
			File dbFolder = new File(dbFolderPath);
			if (!generateHapDatabase(dbFolder, folders[0])) {
				return false;
			}
			updateProgressBar(0.95, 1.0);

			// Step 8: Final checks (100%)
			Platform.runLater(() -> {
				console.appendText("\n✓ All validations passed!\n");
				console.appendText("Database is ready to be added.\n");
			});
			updateProgressBar(1.0, 1.0);

			return true;
		}

		private boolean validateBasicFields() {
			Platform.runLater(() -> console.appendText("\n[Step 1/8] Validating basic fields...\n"));

			if (catalogueId.isEmpty()) {
				Platform.runLater(() -> console.appendText("✗ ERROR: Catalogue ID is required.\n"));
				return false;
			}
			if (identifier.isEmpty()) {
				Platform.runLater(() -> console.appendText("✗ ERROR: Identifier is required.\n"));
				return false;
			}
			if (dbFolderPath.isEmpty()) {
				Platform.runLater(() -> console.appendText("✗ ERROR: Database folder is required.\n"));
				return false;
			}
			if (localVersionTextField == null || localVersionTextField.getText().trim().isEmpty()) {
				Platform.runLater(() -> console.appendText("✗ ERROR: Version name is required.\n"));
				return false;
			}
			if (dbName.isEmpty()) {
				Platform.runLater(() -> console.appendText("✗ ERROR: DB name is required.\n"));
				return false;
			}
			if (funcName.isEmpty()) {
				Platform.runLater(() -> console.appendText("✗ ERROR: Function folder name is required.\n"));
				return false;
			}
			if (taxaName.isEmpty()) {
				Platform.runLater(() -> console.appendText("✗ ERROR: Taxonomy file name is required.\n"));
				return false;
			}

			Platform.runLater(() -> console.appendText("✓ Basic fields validated.\n"));
			return true;
		}

		private boolean validateNoDuplicate() {
			Platform.runLater(() -> console.appendText("\n[Step 1.5/8] Checking for duplicate database...\n"));

			String versionName = localVersionTextField.getText().trim();

			// Check existing local databases in magDbItems array
			if (magDbItems != null && magDbItems[1] != null) {
				for (MagDbItem existingItem : magDbItems[1]) {
					if (existingItem.getCatalogueID().equals(catalogueId)) {
						// Same catalogue ID found, check version
						if (existingItem.getVersionList() != null && !existingItem.getVersionList().isEmpty()) {
							for (MagDbItem.MagVersion existingVersion : existingItem.getVersionList()) {
								if (existingVersion.getCatalogueVersion().equals(versionName)) {
									Platform.runLater(() -> console.appendText("✗ ERROR: A database with Catalogue ID '"
											+ catalogueId + "' and version '" + versionName + "' already exists.\n"));
									Platform.runLater(() -> console
											.appendText("  Please use a different version name or catalogue ID.\n"));
									return false;
								}
							}
						}
					}
				}
			}

			// Also check the table (in case it has items not yet in magDbItems)
			if (localTable != null) {
				for (MagDbItem existingItem : localTable.getItems()) {
					if (existingItem.getCatalogueID().equals(catalogueId)) {
						// Same catalogue ID found, check version
						if (existingItem.getVersionList() != null && !existingItem.getVersionList().isEmpty()) {
							for (MagDbItem.MagVersion existingVersion : existingItem.getVersionList()) {
								if (existingVersion.getCatalogueVersion().equals(versionName)) {
									Platform.runLater(() -> console.appendText("✗ ERROR: A database with Catalogue ID '"
											+ catalogueId + "' and version '" + versionName + "' already exists.\n"));
									Platform.runLater(() -> console
											.appendText("  Please use a different version name or catalogue ID.\n"));
									return false;
								}
							}
						}
					}
				}
			}

			Platform.runLater(() -> console.appendText(
					"✓ No duplicate found. Catalogue ID: " + catalogueId + ", Version: " + versionName + "\n"));
			return true;
		}

		private File[] validateFolderStructure() {
			Platform.runLater(() -> console.appendText("\n[Step 2/8] Validating folder structure...\n"));

			File dbFolder = new File(dbFolderPath);
			if (!dbFolder.exists() || !dbFolder.isDirectory()) {
				Platform.runLater(
						() -> console.appendText("✗ ERROR: Database folder does not exist or is not a directory.\n"));
				return null;
			}
			Platform.runLater(
					() -> console.appendText("✓ Database folder found: " + dbFolder.getAbsolutePath() + "\n"));

			if (dbFolder.getParentFile() == null) {
				Platform.runLater(
						() -> console.appendText("✗ ERROR: Database folder must not be the filesystem root.\n"));
				return null;
			}

			File dbDir = new File(dbFolder, dbName);
			if (!dbDir.exists() || !dbDir.isDirectory()) {
				Platform.runLater(
						() -> console.appendText("✗ ERROR: DB folder not found: " + dbDir.getAbsolutePath() + "\n"));
				return null;
			}
			Platform.runLater(() -> console.appendText("✓ DB folder found: " + dbDir.getAbsolutePath() + "\n"));

			File funcDir = new File(dbFolder, funcName);
			if (!funcDir.exists() || !funcDir.isDirectory()) {
				Platform.runLater(() -> console
						.appendText("✗ ERROR: Function folder not found: " + funcDir.getAbsolutePath() + "\n"));
				return null;
			}
			Platform.runLater(() -> console.appendText("✓ Function folder found: " + funcDir.getAbsolutePath() + "\n"));

			File taxaFile = new File(dbFolder, taxaName);
			if (!taxaFile.exists() || !taxaFile.isFile()) {
				Platform.runLater(() -> console
						.appendText("✗ ERROR: Taxonomy file not found: " + taxaFile.getAbsolutePath() + "\n"));
				return null;
			}
			Platform.runLater(() -> console.appendText("✓ Taxonomy file found: " + taxaFile.getAbsolutePath() + "\n"));

			return new File[] { dbDir, funcDir };
		}

		private boolean validateFileCounts(File originalDbDir, File eggNogDir) {
			Platform.runLater(() -> console.appendText("\n[Step 3/8] Validating file counts...\n"));

			File[] originalDbFiles = originalDbDir.listFiles((dir, name) -> {
				String lower = name.toLowerCase();
				return lower.endsWith(".faa") || lower.endsWith(".fasta");
			});

			File[] eggNogFiles = eggNogDir.listFiles((dir, name) -> name.endsWith(".tsv") || name.endsWith(".txt"));

			if (originalDbFiles == null) {
				Platform.runLater(() -> console.appendText("✗ ERROR: Cannot read original_db directory.\n"));
				return false;
			}
			if (eggNogFiles == null) {
				Platform.runLater(() -> console.appendText("✗ ERROR: Cannot read eggNOG directory.\n"));
				return false;
			}

			Platform.runLater(() -> console.appendText("Found " + originalDbFiles.length + " files in original_db\n"));
			Platform.runLater(() -> console.appendText("Found " + eggNogFiles.length + " files in eggNOG\n"));

			if (originalDbFiles.length != eggNogFiles.length) {
				Platform.runLater(() -> console.appendText("✗ ERROR: File count mismatch! original_db has "
						+ originalDbFiles.length + " files, but eggNOG has " + eggNogFiles.length + " files.\n"));
				return false;
			}

			if (originalDbFiles.length == 0) {
				Platform.runLater(() -> console.appendText("✗ ERROR: No database files found in original_db.\n"));
				return false;
			}

			Platform.runLater(
					() -> console.appendText("✓ File counts match: " + originalDbFiles.length + " pairs found.\n"));

			// Store genome count for later use
			genomeCount = originalDbFiles.length;

			return true;
		}

		private Map<String, File[]> validateFilePairing(File originalDbDir, File eggNogDir) {
			Platform.runLater(() -> console.appendText("\n[Step 4/8] Validating file pairing and naming...\n"));

			File[] originalDbFiles = originalDbDir.listFiles((dir, name) -> {
				String lower = name.toLowerCase();
				return lower.endsWith(".faa") || lower.endsWith(".fasta");
			});

			File[] eggNogFiles = eggNogDir.listFiles((dir, name) -> name.endsWith(".tsv") || name.endsWith(".txt"));

			Map<String, File[]> pairedFiles = new HashMap<>();
			Set<String> processedGenomeIds = new HashSet<>();

			// Process original_db files
			for (File dbFile : originalDbFiles) {
				String fileName = dbFile.getName();
				String genomeId;

				if (fileName.toLowerCase().endsWith(".faa")) {
					genomeId = fileName.substring(0, fileName.length() - 4);
				} else if (fileName.toLowerCase().endsWith(".fasta")) {
					genomeId = fileName.substring(0, fileName.length() - 6);
				} else {
					continue;
				}

				if (genomeId.isEmpty()) {
					Platform.runLater(() -> console.appendText("✗ ERROR: Invalid filename format: " + fileName + "\n"));
					return null;
				}

				// Find matching eggNOG file
				File matchingEggNogFile = null;
				for (File eggNogFile : eggNogFiles) {
					String eggNogName = eggNogFile.getName();
					if (eggNogName.startsWith(genomeId)) {
						matchingEggNogFile = eggNogFile;
						break;
					}
				}

				if (matchingEggNogFile == null) {
					Platform.runLater(() -> console
							.appendText("✗ ERROR: No matching eggNOG file found for genomeID: " + genomeId + "\n"));
					return null;
				}

				if (processedGenomeIds.contains(genomeId)) {
					Platform.runLater(
							() -> console.appendText("✗ ERROR: Duplicate genomeID found: " + genomeId + "\n"));
					return null;
				}

				pairedFiles.put(genomeId, new File[] { dbFile, matchingEggNogFile });
				processedGenomeIds.add(genomeId);
			}

			Platform.runLater(
					() -> console.appendText("✓ All " + pairedFiles.size() + " files are properly paired.\n"));
			return pairedFiles;
		}

		private boolean validateFastaContent(Map<String, File[]> pairedFiles) {
			Platform.runLater(() -> console.appendText("\n[Step 5/8] Validating FASTA content...\n"));

			int totalGenomes = pairedFiles.size();
			final int[] processed = { 0 };
			final int[] errors = { 0 };

			for (Map.Entry<String, File[]> entry : pairedFiles.entrySet()) {
				String genomeId = entry.getKey();
				File fastaFile = entry.getValue()[0];

				processed[0]++;
				final int currentProcessed = processed[0];
				double progress = 0.5 + (0.3 * currentProcessed / totalGenomes);
				updateProgressBar(progress, 1.0);

				Platform.runLater(() -> console.appendText(
						"Validating " + genomeId + " (" + currentProcessed + "/" + totalGenomes + ")...\n"));

				try (BufferedReader reader = new BufferedReader(new FileReader(fastaFile))) {
					String line;
					final int[] proteinCount = { 0 };
					boolean hasError = false;

					while ((line = reader.readLine()) != null) {
						if (line.startsWith(">")) {
							proteinCount[0]++;

							// Check format: >genomeID_proteinID protein name
							String header = line.substring(1).trim();

							if (!header.startsWith(genomeId + "_")) {
								Platform.runLater(() -> console.appendText("✗ ERROR in " + genomeId
										+ ": Protein header does not start with genomeID: " + header + "\n"));
								hasError = true;
								errors[0]++;
								break;
							}

							// Check for proteinID and protein name
							String afterGenomeId = header.substring(genomeId.length() + 1);
							int spaceIndex = afterGenomeId.indexOf(' ');

							if (spaceIndex <= 0) {
								Platform.runLater(() -> console.appendText("✗ ERROR in " + genomeId
										+ ": Missing proteinID or protein name: " + header + "\n"));
								hasError = true;
								errors[0]++;
								break;
							}

							String proteinId = afterGenomeId.substring(0, spaceIndex);
							String proteinName = afterGenomeId.substring(spaceIndex + 1).trim();

							if (proteinId.isEmpty()) {
								Platform.runLater(() -> console
										.appendText("✗ ERROR in " + genomeId + ": Empty proteinID: " + header + "\n"));
								hasError = true;
								errors[0]++;
								break;
							}

							if (proteinName.isEmpty()) {
								Platform.runLater(() -> console.appendText("✗ ERROR in " + genomeId
										+ ": Missing protein name (mandatory): " + header + "\n"));
								hasError = true;
								errors[0]++;
								break;
							}
						}
					}

					if (!hasError) {
						final int finalProteinCount = proteinCount[0];
						Platform.runLater(() -> console
								.appendText("  ✓ " + genomeId + ": " + finalProteinCount + " proteins validated\n"));
					}

					if (errors[0] > 10) {
						Platform.runLater(() -> console.appendText("✗ Too many errors. Stopping validation.\n"));
						return false;
					}

				} catch (IOException e) {
					Platform.runLater(() -> console
							.appendText("✗ ERROR reading " + fastaFile.getName() + ": " + e.getMessage() + "\n"));
					return false;
				}
			}

			if (errors[0] > 0) {
				final int finalErrors = errors[0];
				Platform.runLater(
						() -> console.appendText("\n✗ Validation failed with " + finalErrors + " error(s).\n"));
				return false;
			}

			Platform.runLater(() -> console.appendText("✓ All FASTA files validated successfully.\n"));
			return true;
		}

		private boolean validateEggNogFormat(Map<String, File[]> pairedFiles, File eggNogDir) {
			Platform.runLater(() -> console.appendText("\n[Step 6/8] Validating eggNOG format...\n"));

			// Basic validation: check that eggNOG files exist and are readable
			// eggnogMapper format typically has tab-separated values with specific columns
			final int[] validated = { 0 };
			int total = pairedFiles.size();

			for (Map.Entry<String, File[]> entry : pairedFiles.entrySet()) {
				File eggNogFile = entry.getValue()[1];

				validated[0]++;
				final int currentValidated = validated[0];
				double progress = 0.8 + (0.15 * currentValidated / total);
				updateProgressBar(progress, 1.0);

				if (!eggNogFile.exists() || !eggNogFile.canRead()) {
					Platform.runLater(() -> console
							.appendText("✗ ERROR: Cannot read eggNOG file: " + eggNogFile.getName() + "\n"));
					return false;
				}

				// Basic format check: should have tab-separated columns
				try (BufferedReader reader = new BufferedReader(new FileReader(eggNogFile))) {
					String firstLine = reader.readLine();
					if (firstLine == null || firstLine.trim().isEmpty()) {
						Platform.runLater(
								() -> console.appendText("✗ ERROR: Empty eggNOG file: " + eggNogFile.getName() + "\n"));
						return false;
					}

					// Check if it's tab-separated (eggnogMapper format)
					String[] columns = firstLine.split("\t");
					if (columns.length < 2) {
						Platform.runLater(
								() -> console.appendText("⚠ WARNING: eggNOG file may not be in correct format: "
										+ eggNogFile.getName() + " (only " + columns.length + " columns found)\n"));
					}

					// Check that file has some content
					int lineCount = 1;
					while (reader.readLine() != null) {
						lineCount++;
						if (lineCount > 10)
							break; // Just check first few lines
					}

					if (lineCount < 2) {
						Platform.runLater(() -> console.appendText(
								"⚠ WARNING: eggNOG file has very few lines: " + eggNogFile.getName() + "\n"));
					}

				} catch (IOException e) {
					Platform.runLater(() -> console.appendText(
							"✗ ERROR reading eggNOG file " + eggNogFile.getName() + ": " + e.getMessage() + "\n"));
					return false;
				}
			}

			final int finalValidated = validated[0];
			Platform.runLater(
					() -> console.appendText("✓ All eggNOG files validated (" + finalValidated + "/" + total + ").\n"));
			return true;
		}

		private boolean generateHapDatabase(File dbFolder, File originalDbDir) {
			Platform.runLater(() -> console.appendText("\n[Step 7/8] Generating HAP database...\n"));

			try {
				// 4-markers: Anchor, HighRes, Abundance, Stability (use extractFourCategories)
				HapExtractor extractor = new HapExtractor();

				File marker2File = new File(dbFolder, "Marker2.fasta");
				File marker4File = new File(dbFolder, "Marker4.fasta");

				Platform.runLater(() -> console.appendText("Extracting Marker2 proteins from original_db...\n"));
				Platform.runLater(() -> console.appendText("Output file: " + marker2File.getAbsolutePath() + "\n"));

				extractor.extractTwoMarkers(originalDbDir, marker2File);

				Platform.runLater(() -> console.appendText("Extracting Marker4 proteins from original_db...\n"));
				Platform.runLater(() -> console.appendText("Output file: " + marker4File.getAbsolutePath() + "\n"));

				extractor.extractFourCategories(originalDbDir, marker4File);

				// Verify the file was created and is not empty
				if (!marker2File.exists()) {
					Platform.runLater(() -> console.appendText(
							"✗ ERROR: Marker2 file was not created: " + marker2File.getAbsolutePath() + "\n"));
					return false;
				}

				// Verify the file was created and is not empty
				if (!marker4File.exists()) {
					Platform.runLater(() -> console.appendText(
							"✗ ERROR: Marker4 file was not created: " + marker4File.getAbsolutePath() + "\n"));
					return false;
				}

				if (marker2File.length() == 0) {
					Platform.runLater(() -> console
							.appendText("✗ ERROR: Marker2 file is empty: " + marker2File.getAbsolutePath() + "\n"));
					marker2File.delete(); // Clean up empty file
					return false;
				}

				if (marker4File.length() == 0) {
					Platform.runLater(() -> console
							.appendText("✗ ERROR: Marker4 file is empty: " + marker4File.getAbsolutePath() + "\n"));
					marker4File.delete(); // Clean up empty file
					return false;
				}

				// Check if file has content (at least some headers)
				final int[] headerCount = new int[2];
				try (BufferedReader reader = new BufferedReader(new FileReader(marker2File))) {
					// Re-read to count headers
					String line;
					while ((line = reader.readLine()) != null) {
						if (line.startsWith(">")) {
							headerCount[0]++;
						}
					}
				}

				if (headerCount[0] == 0) {
					Platform.runLater(
							() -> console.appendText("✗ ERROR: Marker2 file contains no protein sequences\n"));
					marker2File.delete(); // Clean up invalid file
					return false;
				}

				try (BufferedReader reader = new BufferedReader(new FileReader(marker4File))) {
					// Re-read to count headers
					String line;
					while ((line = reader.readLine()) != null) {
						if (line.startsWith(">")) {
							headerCount[1]++;
						}
					}
				}

				if (headerCount[1] == 0) {
					Platform.runLater(
							() -> console.appendText("✗ ERROR: Marker4 file contains no protein sequences\n"));
					marker4File.delete(); // Clean up invalid file
					return false;
				}

				Platform.runLater(() -> console.appendText("✓ Marker2 and Marker4 database generated successfully\n"));
				Platform.runLater(() -> console.appendText("  File size: " + marker2File.length() + " bytes\n"
						+ "  File size: \" + marker4File.length() + \" bytes"));
				Platform.runLater(() -> console
						.appendText("  Found at least " + headerCount[0] + " protein sequences in Marker2\n"
								+ "  Found at least \" + headerCount[1] + \" protein sequences in Marker4"));

				return true;

			} catch (IOException e) {
				Platform.runLater(
						() -> console.appendText("✗ ERROR: Failed to generate HAP database: " + e.getMessage() + "\n"));
				e.printStackTrace();
				return false;
			} catch (Exception e) {
				Platform.runLater(() -> console
						.appendText("✗ ERROR: Unexpected error during HAP generation: " + e.getMessage() + "\n"));
				e.printStackTrace();
				return false;
			}
		}

		private void updateProgressBar(double workDone, double max) {
			Platform.runLater(() -> {
				progressBar.setProgress(workDone / max);
			});
		}
	}

	private MagDbItem buildLocalMagDbItem(int genomeCount) {
		String catalogueId = localNameTextField.getText().trim();
		String identifier = identifierTextField.getText().trim();
		File dbFolder = new File(localRepoTextField.getText().trim());
		// For local databases, the repository folder is the database folder itself
		File repositoryFolder = dbFolder;
		String versionName = localVersionTextField.getText().trim();

		String dbName = localDbTextField.getText().trim();
		String funcName = localFuncTextField.getText().trim();
		String taxaName = localTaxaTextField.getText().trim();

		MagDbItem newItem = new MagDbItem(catalogueId, identifier, repositoryFolder);
		newItem.addVersion(versionName, genomeCount, new Date(), true, dbName, funcName, taxaName);
		newItem.setUsedVersion(newItem.getVersionList().get(0));

		// For local databases, the currentFile should be the database folder itself
		// (not repositoryFolder/versionName as setUsedVersion would create)
		newItem.setCurrentFile(dbFolder);

		return newItem;
	}

	private void addLocalMagDbItem(MagDbItem newItem) {
		if (magDbItems == null) {
			magDbItems = new MagDbItem[2][];
		}
		if (magDbItems[1] == null) {
			magDbItems[1] = new MagDbItem[] { newItem };
		} else {
			MagDbItem[] newArray = new MagDbItem[magDbItems[1].length + 1];
			System.arraycopy(magDbItems[1], 0, newArray, 0, magDbItems[1].length);
			newArray[magDbItems[1].length] = newItem;
			magDbItems[1] = newArray;
		}

		// Save to file
		saveMagDbConfig();
	}

	private void applyLocalMagDbItem(MagDbItem newItem) {
		if (newItem == null) {
			return;
		}

		// Add to magDbItems array first
		if (magDbItems == null) {
			magDbItems = new MagDbItem[2][];
		}
		if (magDbItems[1] == null) {
			magDbItems[1] = new MagDbItem[] { newItem };
		} else {
			MagDbItem[] newArray = new MagDbItem[magDbItems[1].length + 1];
			System.arraycopy(magDbItems[1], 0, newArray, 0, magDbItems[1].length);
			newArray[magDbItems[1].length] = newItem;
			magDbItems[1] = newArray;
		}

		// Update table
		loadLocalTableData();

		// Save to file (after table is updated)
		saveMagDbConfig();
	}

	/**
	 * Save the current MAG database configuration to file
	 */
	private void saveMagDbConfig() {
		try {
			// Update local database array from table
			if (localTable != null && !localTable.getItems().isEmpty()) {
				magDbItems[1] = localTable.getItems().toArray(new MagDbItem[0]);
			}

			// Get the config file path (same as used in MetaPilotApp)
			File dbConfigFile = new File(new File("").getAbsolutePath(), MagDbConfigIO.versionFile);

			// Ensure magDbItems array is properly initialized
			if (magDbItems == null) {
				magDbItems = new MagDbItem[2][];
			}
			if (magDbItems[0] == null) {
				magDbItems[0] = new MagDbItem[0];
			}
			if (magDbItems[1] == null) {
				magDbItems[1] = new MagDbItem[0];
			}

			// Export to file
			MagDbConfigIO.export(magDbItems, dbConfigFile);

		} catch (Exception e) {
			showAlert(Alert.AlertType.ERROR, "Save Failed",
					"Failed to save MAG database configuration: " + e.getMessage());
		}
	}

	/**
	 * Helper method to show alerts
	 */
	private void showAlert(Alert.AlertType type, String title, String message) {
		Alert alert = new Alert(type);
		alert.setTitle(title);
		alert.setHeaderText(null);
		alert.setContentText(message);
		alert.showAndWait();
	}

	/**
	 * Get the current database configuration
	 */
	public MagDbItem[][] getMagDbItems() {
		// Update local database array from table
		if (localTable != null && !localTable.getItems().isEmpty()) {
			magDbItems[1] = localTable.getItems().toArray(new MagDbItem[0]);
		}

		return magDbItems;
	}

}
