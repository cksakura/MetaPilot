package uk.ac.quadram.metapilot.ui;

import java.io.File;
import java.nio.file.Path;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.List;
import java.util.Optional;

import javafx.beans.property.ReadOnlyStringWrapper;
import javafx.collections.FXCollections;
import javafx.collections.ObservableList;
import javafx.geometry.Insets;
import javafx.geometry.Pos;
import javafx.scene.control.Alert;
import javafx.scene.control.Button;
import javafx.scene.control.ButtonType;
import javafx.scene.control.CheckBox;
import javafx.scene.control.ComboBox;
import javafx.scene.control.Label;
import javafx.scene.control.RadioButton;
import javafx.scene.control.ScrollPane;
import javafx.scene.control.TableColumn;
import javafx.scene.control.TableView;
import javafx.scene.control.TextArea;
import javafx.scene.control.TextField;
import javafx.scene.control.ToggleGroup;
import javafx.scene.layout.BorderPane;
import javafx.scene.layout.GridPane;
import javafx.scene.layout.HBox;
import javafx.scene.layout.Priority;
import javafx.scene.layout.VBox;
import javafx.stage.DirectoryChooser;

import uk.ac.quadram.metapilot.core.clean.ResultFolderCleaner;
import uk.ac.quadram.metapilot.core.clean.ResultFolderCleaner.CleanMode;
import uk.ac.quadram.metapilot.core.clean.ResultFolderCleaner.CleanPreview;
import uk.ac.quadram.metapilot.core.clean.ResultFolderCleaner.CleanResult;
import uk.ac.quadram.metapilot.core.clean.ResultFolderCleaner.DetectionResult;
import uk.ac.quadram.metapilot.core.clean.ResultFolderCleaner.PreviewEntry;
import uk.ac.quadram.metapilot.workflow.WorkflowType;

/**
 * JavaFX panel for previewing and cleaning MetaPilot result folders.
 */
public class ResultCleanerPanel extends BorderPane {

	private static final DecimalFormat SIZE_FORMAT = new DecimalFormat("#,##0.0");

	private final ResultFolderCleaner cleaner = new ResultFolderCleaner();
	private final TextField resultFolderField = new TextField();
	private final Label detectionLabel = new Label("Select a result folder to detect workflow.");
	private final CheckBox manualWorkflowCheckBox = new CheckBox("Choose workflow manually");
	private final ComboBox<WorkflowChoice> workflowComboBox = new ComboBox<>();
	private final ComboBox<String> diaModeComboBox = new ComboBox<>();
	private final ToggleGroup cleanModeGroup = new ToggleGroup();
	private final RadioButton basicCleanRadio = new RadioButton("Basic clean");
	private final RadioButton deepCleanRadio = new RadioButton("Deep clean");
	private final CheckBox removeConvertedSpectraCheckBox = new CheckBox("Remove converted spectra files");
	private final Button previewButton = new Button("Preview Clean");
	private final Button confirmButton = new Button("Confirm Clean");
	private final TextArea summaryArea = new TextArea();
	private final TableView<PreviewRow> previewTable = new TableView<>();
	private final ObservableList<PreviewRow> previewRows = FXCollections.observableArrayList();

	private DetectionResult detectionResult;
	private CleanPreview currentPreview;

	public ResultCleanerPanel() {
		setPadding(new Insets(20));
		getStyleClass().add("section-panel");
		setCenter(createContent());
		initializeState();
	}

	private ScrollPane createContent() {
		VBox content = new VBox(18);
		content.setFillWidth(true);

		Label title = new Label("Result Folder Cleaner");
		title.getStyleClass().add("section-title");
		title.setStyle("-fx-font-size: 20px;");

		Label intro = wrapLabel(
				"Clean completed MetaPilot result folders with a preview-first workflow. Basic clean keeps final outputs plus rerun/resume checkpoints. Deep clean keeps final deliverables, with converted spectra controlled by a separate option.");

		content.getChildren().addAll(title, intro, createFolderSection(), createWorkflowSection(), createModeSection(),
				createPreviewSection());

		ScrollPane scrollPane = new ScrollPane(content);
		scrollPane.setFitToWidth(true);
		scrollPane.setPadding(new Insets(0));
		return scrollPane;
	}

	private VBox createFolderSection() {
		VBox section = sectionBox("1. Result folder");
		Label help = wrapLabel("Select the root result folder produced by one MetaPilot workflow run.");

		Button browseButton = new Button("Browse...");
		browseButton.setOnAction(e -> browseForFolder());

		Button detectButton = new Button("Detect Workflow");
		detectButton.setOnAction(e -> detectWorkflow());

		resultFolderField.setPromptText("Path to MetaPilot result folder");
		resultFolderField.textProperty().addListener((obs, oldValue, newValue) -> {
			clearPreview();
			detectionResult = null;
			detectionLabel.setText("Workflow detection has not been run for this folder.");
			workflowComboBox.getSelectionModel().selectFirst();
			workflowComboBox.setDisable(!manualWorkflowCheckBox.isSelected());
			updateActionState();
		});
		HBox.setHgrow(resultFolderField, Priority.ALWAYS);

		HBox row = new HBox(10, resultFolderField, browseButton, detectButton);
		row.setAlignment(Pos.CENTER_LEFT);
		section.getChildren().addAll(help, row);
		return section;
	}

	private VBox createWorkflowSection() {
		VBox section = sectionBox("2. Workflow");
		Label help = wrapLabel(
				"The cleaner first reads saved MetaPilot parameter files to identify the workflow. If the result folder is old, copied, or ambiguous, enable manual selection.");

		detectionLabel.setWrapText(true);
		detectionLabel.setStyle("-fx-font-weight: bold;");

		workflowComboBox.setItems(FXCollections.observableArrayList(createWorkflowChoices()));
		workflowComboBox.setMaxWidth(Double.MAX_VALUE);
		workflowComboBox.valueProperty().addListener((obs, oldValue, newValue) -> {
			updateDiaModeVisibility();
			clearPreview();
			updateActionState();
		});

		diaModeComboBox.setItems(FXCollections.observableArrayList("Auto / not specified", "Direct", "Core / Focus",
				"Comprehensive", "Deep together", "Deep per file"));
		diaModeComboBox.getSelectionModel().selectFirst();
		diaModeComboBox.valueProperty().addListener((obs, oldValue, newValue) -> clearPreview());

		manualWorkflowCheckBox.selectedProperty().addListener((obs, oldValue, selected) -> {
			workflowComboBox.setDisable(!selected && detectionResult != null && detectionResult.hasWorkflow());
			updateActionState();
		});

		GridPane grid = new GridPane();
		grid.setHgap(12);
		grid.setVgap(10);
		grid.add(new Label("Detected:"), 0, 0);
		grid.add(detectionLabel, 1, 0);
		grid.add(manualWorkflowCheckBox, 1, 1);
		grid.add(new Label("Workflow:"), 0, 2);
		grid.add(workflowComboBox, 1, 2);
		grid.add(new Label("DIA mode:"), 0, 3);
		grid.add(diaModeComboBox, 1, 3);
		GridPane.setHgrow(workflowComboBox, Priority.ALWAYS);
		GridPane.setHgrow(detectionLabel, Priority.ALWAYS);

		section.getChildren().addAll(help, grid);
		return section;
	}

	private VBox createModeSection() {
		VBox section = sectionBox("3. Clean mode");

		basicCleanRadio.setToggleGroup(cleanModeGroup);
		deepCleanRadio.setToggleGroup(cleanModeGroup);
		basicCleanRadio.setSelected(true);
		cleanModeGroup.selectedToggleProperty().addListener((obs, oldValue, newValue) -> {
			clearPreview();
			updateActionState();
		});

		VBox basicBox = modeBox(basicCleanRadio, CleanMode.BASIC.getDescription());
		VBox deepBox = modeBox(deepCleanRadio, CleanMode.DEEP.getDescription());

		HBox modes = new HBox(16, basicBox, deepBox);
		modes.setFillHeight(true);
		HBox.setHgrow(basicBox, Priority.ALWAYS);
		HBox.setHgrow(deepBox, Priority.ALWAYS);

		Label finalRules = wrapLabel(
				"Final deliverables include top-level peptide, protein, summary and report files; taxonomy or taxonomy_analysis folders; function or functional_annotation folders; report folders; and DIA combined outputs when the selected workflow is DIA.");

		VBox spectraBox = createConvertedSpectraSection();

		section.getChildren().addAll(modes, finalRules, spectraBox);
		return section;
	}

	private VBox createConvertedSpectraSection() {
		VBox spectraBox = new VBox(8);
		spectraBox.setPadding(new Insets(12));
		spectraBox.getStyleClass().add("workflow-box");

		Label title = new Label("Converted spectra");
		title.setStyle("-fx-font-weight: bold;");
		Label description = wrapLabel(
				"Converted spectra are large intermediate files produced early in workflows such as pFind. Keep them when you want Basic clean to preserve rerun/resume behavior. Removing them saves space, but reruns will start from spectra conversion.");

		removeConvertedSpectraCheckBox.setSelected(false);
		removeConvertedSpectraCheckBox.setOnAction(e -> {
			clearPreview();
			updateActionState();
		});

		spectraBox.getChildren().addAll(title, description, removeConvertedSpectraCheckBox);
		return spectraBox;
	}

	private VBox createPreviewSection() {
		VBox section = sectionBox("4. Preview and confirm");
		Label help = wrapLabel(
				"Preview lists every retained and deleted path. The confirm button is enabled only after the current folder, workflow and clean mode have been previewed.");

		previewButton.setOnAction(e -> buildPreview());
		confirmButton.setOnAction(e -> confirmClean());
		confirmButton.setDisable(true);

		HBox actions = new HBox(10, previewButton, confirmButton);
		actions.setAlignment(Pos.CENTER_LEFT);

		summaryArea.setEditable(false);
		summaryArea.setWrapText(true);
		summaryArea.setPrefRowCount(5);
		summaryArea.setText("No preview yet.");

		configurePreviewTable();

		section.getChildren().addAll(help, actions, summaryArea, previewTable);
		VBox.setVgrow(previewTable, Priority.ALWAYS);
		return section;
	}

	private void configurePreviewTable() {
		TableColumn<PreviewRow, String> actionColumn = new TableColumn<>("Action");
		actionColumn.setCellValueFactory(data -> new ReadOnlyStringWrapper(data.getValue().action()));
		actionColumn.setPrefWidth(80);

		TableColumn<PreviewRow, String> typeColumn = new TableColumn<>("Type");
		typeColumn.setCellValueFactory(data -> new ReadOnlyStringWrapper(data.getValue().type()));
		typeColumn.setPrefWidth(80);

		TableColumn<PreviewRow, String> pathColumn = new TableColumn<>("Path");
		pathColumn.setCellValueFactory(data -> new ReadOnlyStringWrapper(data.getValue().path()));
		pathColumn.setPrefWidth(520);

		TableColumn<PreviewRow, String> sizeColumn = new TableColumn<>("Size");
		sizeColumn.setCellValueFactory(data -> new ReadOnlyStringWrapper(data.getValue().size()));
		sizeColumn.setPrefWidth(100);

		TableColumn<PreviewRow, String> reasonColumn = new TableColumn<>("Reason");
		reasonColumn.setCellValueFactory(data -> new ReadOnlyStringWrapper(data.getValue().reason()));
		reasonColumn.setPrefWidth(360);

		previewTable.getColumns().addAll(actionColumn, typeColumn, pathColumn, sizeColumn, reasonColumn);
		previewTable.setItems(previewRows);
		previewTable.setPrefHeight(380);
		previewTable.setColumnResizePolicy(TableView.CONSTRAINED_RESIZE_POLICY_FLEX_LAST_COLUMN);
	}

	private VBox sectionBox(String titleText) {
		VBox section = new VBox(10);
		section.getStyleClass().add("section-panel");
		Label title = new Label(titleText);
		title.getStyleClass().add("section-title");
		section.getChildren().add(title);
		return section;
	}

	private VBox modeBox(RadioButton radioButton, String description) {
		Label descriptionLabel = wrapLabel(description);
		VBox box = new VBox(8, radioButton, descriptionLabel);
		box.setPadding(new Insets(12));
		box.getStyleClass().add("workflow-box");
		return box;
	}

	private Label wrapLabel(String text) {
		Label label = new Label(text);
		label.setWrapText(true);
		return label;
	}

	private void initializeState() {
		workflowComboBox.getSelectionModel().selectFirst();
		workflowComboBox.setDisable(true);
		diaModeComboBox.setDisable(true);
		updateActionState();
	}

	private void browseForFolder() {
		DirectoryChooser chooser = new DirectoryChooser();
		chooser.setTitle("Select MetaPilot result folder");
		String currentPath = resultFolderField.getText();
		if (currentPath != null && !currentPath.isBlank()) {
			File current = new File(currentPath);
			if (current.isDirectory()) {
				chooser.setInitialDirectory(current);
			}
		}
		File selected = chooser.showDialog(getScene() != null ? getScene().getWindow() : null);
		if (selected != null) {
			resultFolderField.setText(selected.getAbsolutePath());
			detectWorkflow();
		}
	}

	private void detectWorkflow() {
		clearPreview();
		Path root = selectedRoot();
		if (root == null) {
			showAlert(Alert.AlertType.WARNING, "Result Folder Required", "Select an existing result folder first.");
			return;
		}

		detectionResult = cleaner.detect(root);
		if (detectionResult.hasWorkflow()) {
			selectWorkflow(detectionResult.getWorkflowType());
			detectionLabel.setText("Detected " + displayWorkflow(detectionResult.getWorkflowType())
					+ (detectionResult.isHighConfidence() ? " from saved parameter metadata." : " from available markers."));
			if (detectionResult.getRefinementMode() != null) {
				diaModeComboBox.getSelectionModel().select(detectionResult.getRefinementMode());
			}
			workflowComboBox.setDisable(!manualWorkflowCheckBox.isSelected());
		} else {
			detectionLabel.setText("Workflow was not detected. Choose workflow manually before previewing.");
			manualWorkflowCheckBox.setSelected(true);
			workflowComboBox.getSelectionModel().selectFirst();
			workflowComboBox.setDisable(false);
		}
		updateDiaModeVisibility();
		updateSummaryForDetection();
		updateActionState();
	}

	private void buildPreview() {
		Path root = selectedRoot();
		WorkflowType workflowType = selectedWorkflow();
		if (root == null) {
			showAlert(Alert.AlertType.WARNING, "Result Folder Required", "Select an existing result folder first.");
			return;
		}
		if (workflowType == null) {
			showAlert(Alert.AlertType.WARNING, "Workflow Required", "Detect or choose the workflow before previewing.");
			return;
		}

		try {
			currentPreview = cleaner.preview(root, workflowType, selectedDiaMode(), selectedCleanMode(),
					removeConvertedSpectraCheckBox.isSelected());
			previewRows.setAll(createPreviewRows(currentPreview));
			updateSummaryForPreview(currentPreview);
			showPreviewWarnings(currentPreview);
			confirmButton.setDisable(false);
		} catch (Exception ex) {
			clearPreview();
			String message = ex.getMessage() != null ? ex.getMessage() : "Preview failed.";
			if (message.startsWith("The selected folder")) {
				showAlert(Alert.AlertType.WARNING, "Root Folder Not Accepted", message);
			} else {
				showAlert(Alert.AlertType.ERROR, "Preview Failed", message);
			}
		}
	}

	private void confirmClean() {
		if (currentPreview == null) {
			return;
		}
		Alert confirm = new Alert(Alert.AlertType.CONFIRMATION);
		confirm.setTitle("Confirm Result Folder Clean");
		confirm.setHeaderText("Delete " + currentPreview.getFilesToDelete() + " files and "
				+ currentPreview.getDirectoriesToDelete() + " folders?");
		String spectraWarning = currentPreview.isRemoveConvertedSpectra()
				? "\n\nConverted spectra are included in the delete preview. Reruns will need to perform spectra conversion again."
				: "";
		confirm.setContentText("This will delete only the paths shown as DELETE in the preview table."
				+ spectraWarning + "\n\nResult folder:\n" + currentPreview.getResultRoot());
		Optional<ButtonType> response = confirm.showAndWait();
		if (response.isEmpty() || response.get() != ButtonType.OK) {
			return;
		}

		try {
			CleanResult result = cleaner.clean(currentPreview);
			showAlert(Alert.AlertType.INFORMATION, "Clean Complete",
					"Deleted " + result.getDeletedFiles() + " files and " + result.getDeletedDirectories()
							+ " folders.\nSpace removed: " + formatBytes(result.getDeletedBytes()) + "\nLog file: "
							+ result.getLogFile());
			clearPreview();
			detectWorkflow();
		} catch (Exception ex) {
			showAlert(Alert.AlertType.ERROR, "Clean Failed", ex.getMessage());
		}
	}

	private List<PreviewRow> createPreviewRows(CleanPreview preview) {
		List<PreviewRow> rows = new ArrayList<>();
		for (PreviewEntry entry : preview.getEntries()) {
			Path relative = preview.getResultRoot().relativize(entry.getPath());
			rows.add(new PreviewRow(entry.isKeep() ? "KEEP" : "DELETE", entry.isDirectory() ? "Folder" : "File",
					relative.toString(), entry.isDirectory() ? "" : formatBytes(entry.getSizeBytes()), entry.getReason()));
		}
		return rows;
	}

	private void clearPreview() {
		currentPreview = null;
		previewRows.clear();
		confirmButton.setDisable(true);
		summaryArea.setText("No preview yet.");
	}

	private void updateSummaryForDetection() {
		if (detectionResult == null) {
			return;
		}
		StringBuilder sb = new StringBuilder();
		if (!detectionResult.getEvidence().isEmpty()) {
			sb.append("Detection evidence:\n");
			for (String evidence : detectionResult.getEvidence()) {
				sb.append("- ").append(evidence).append("\n");
			}
		}
		if (!detectionResult.getWarnings().isEmpty()) {
			if (sb.length() > 0) {
				sb.append("\n");
			}
			sb.append("Warnings:\n");
			for (String warning : detectionResult.getWarnings()) {
				sb.append("- ").append(warning).append("\n");
			}
		}
		if (sb.length() == 0) {
			sb.append("No preview yet.");
		}
		summaryArea.setText(sb.toString());
	}

	private void updateSummaryForPreview(CleanPreview preview) {
		StringBuilder sb = new StringBuilder();
		sb.append("Preview for ").append(displayWorkflow(preview.getWorkflowType())).append(" / ")
				.append(preview.getMode().getDisplayName()).append("\n");
		sb.append("Files to delete: ").append(preview.getFilesToDelete()).append("\n");
		sb.append("Folders to delete: ").append(preview.getDirectoriesToDelete()).append("\n");
		sb.append("Estimated space saved: ").append(formatBytes(preview.getBytesToDelete())).append("\n");
		sb.append("Converted spectra: ")
				.append(preview.isRemoveConvertedSpectra() ? "selected for removal" : "kept for rerun/resume")
				.append("\n");
		if (!preview.getWarnings().isEmpty()) {
			sb.append("\nWarnings:\n");
			for (String warning : preview.getWarnings()) {
				sb.append("- ").append(warning).append("\n");
			}
		}
		summaryArea.setText(sb.toString());
	}

	private void showPreviewWarnings(CleanPreview preview) {
		if (preview == null || preview.getWarnings().isEmpty()) {
			return;
		}
		StringBuilder sb = new StringBuilder();
		for (String warning : preview.getWarnings()) {
			if (sb.length() > 0) {
				sb.append("\n\n");
			}
			sb.append(warning);
		}
		showAlert(Alert.AlertType.WARNING, "Review Preview Warnings", sb.toString());
	}

	private void updateActionState() {
		previewButton.setDisable(selectedRoot() == null || selectedWorkflow() == null || selectedCleanMode() == null);
	}

	private void updateDiaModeVisibility() {
		WorkflowType workflowType = selectedWorkflow();
		boolean diaMagMode = workflowType == WorkflowType.DIA_MARKER || workflowType == WorkflowType.DIA_SEED;
		diaModeComboBox.setDisable(!diaMagMode);
		if (!diaMagMode) {
			diaModeComboBox.getSelectionModel().selectFirst();
		}
	}

	private Path selectedRoot() {
		String path = resultFolderField.getText();
		if (path == null || path.isBlank()) {
			return null;
		}
		File file = new File(path.trim());
		if (!file.isDirectory()) {
			return null;
		}
		return file.toPath();
	}

	private CleanMode selectedCleanMode() {
		if (cleanModeGroup.getSelectedToggle() == basicCleanRadio) {
			return CleanMode.BASIC;
		}
		if (cleanModeGroup.getSelectedToggle() == deepCleanRadio) {
			return CleanMode.DEEP;
		}
		return null;
	}

	private WorkflowType selectedWorkflow() {
		WorkflowChoice choice = workflowComboBox.getSelectionModel().getSelectedItem();
		return choice != null ? choice.workflowType() : null;
	}

	private String selectedDiaMode() {
		if (diaModeComboBox.isDisabled()) {
			return null;
		}
		String mode = diaModeComboBox.getSelectionModel().getSelectedItem();
		if (mode == null || mode.startsWith("Auto")) {
			return detectionResult != null ? detectionResult.getRefinementMode() : null;
		}
		return mode;
	}

	private void selectWorkflow(WorkflowType workflowType) {
		for (WorkflowChoice choice : workflowComboBox.getItems()) {
			if (choice.workflowType() == workflowType) {
				workflowComboBox.getSelectionModel().select(choice);
				return;
			}
		}
	}

	private List<WorkflowChoice> createWorkflowChoices() {
		return List.of(
				new WorkflowChoice("Select workflow...", null),
				new WorkflowChoice("DDA / IGC Database / MaxQuant", WorkflowType.DDA_IGC_MAXQUANT),
				new WorkflowChoice("DDA / IGC Database / pFind", WorkflowType.DDA_IGC_PFIND),
				new WorkflowChoice("DDA / IGC Database / FragPipe", WorkflowType.DDA_IGC_FRAGPIPE),
				new WorkflowChoice("DDA / IGC Database / Sage", WorkflowType.DDA_IGC_SAGE),
				new WorkflowChoice("DDA / MAG Database / pFind", WorkflowType.DDA_MAG_PFIND),
				new WorkflowChoice("DDA / MAG Database / FragPipe", WorkflowType.DDA_MAG_FRAGPIPE),
				new WorkflowChoice("DDA / MAG Database / Sage", WorkflowType.DDA_MAG_SAGE),
				new WorkflowChoice("DIA / Marker", WorkflowType.DIA_MARKER),
				new WorkflowChoice("DIA / Seeded", WorkflowType.DIA_SEED),
				new WorkflowChoice("DIA / Target", WorkflowType.DIA_TARGET),
				new WorkflowChoice("DIA / Real-time", WorkflowType.DIA_REALTIME));
	}

	private String displayWorkflow(WorkflowType workflowType) {
		if (workflowType == null) {
			return "Unknown workflow";
		}
		for (WorkflowChoice choice : workflowComboBox.getItems()) {
			if (choice.workflowType() == workflowType) {
				return choice.label();
			}
		}
		return workflowType.getAnalysisType();
	}

	private String formatBytes(long bytes) {
		if (bytes < 1024L) {
			return bytes + " B";
		}
		double kb = bytes / 1024.0;
		if (kb < 1024.0) {
			return SIZE_FORMAT.format(kb) + " KB";
		}
		double mb = kb / 1024.0;
		if (mb < 1024.0) {
			return SIZE_FORMAT.format(mb) + " MB";
		}
		return SIZE_FORMAT.format(mb / 1024.0) + " GB";
	}

	private void showAlert(Alert.AlertType type, String title, String message) {
		Alert alert = new Alert(type);
		alert.setTitle(title);
		alert.setHeaderText(null);
		alert.setContentText(message);
		alert.showAndWait();
	}

	private record WorkflowChoice(String label, WorkflowType workflowType) {
		@Override
		public String toString() {
			return label;
		}
	}

	private record PreviewRow(String action, String type, String path, String size, String reason) {
	}
}
