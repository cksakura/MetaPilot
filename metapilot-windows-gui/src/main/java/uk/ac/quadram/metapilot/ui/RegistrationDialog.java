package uk.ac.quadram.metapilot.ui;

import javafx.application.Platform;
import javafx.geometry.Insets;
import javafx.geometry.Pos;
import javafx.scene.control.*;
import javafx.scene.layout.*;
import javafx.stage.Stage;
import javafx.stage.WindowEvent;
import javafx.scene.Node;

import java.util.Arrays;
import java.util.List;
import java.util.Locale;
import java.awt.Desktop;
import java.net.URI;

import uk.ac.quadram.metapilot.registration.RegistrationManager;
import uk.ac.quadram.metapilot.registration.UserRegistration;

/**
 * Dialog for user registration and activation
 */
public class RegistrationDialog extends Dialog<UserRegistration> {
    
    private TextField nameField;
    private TextField emailField;
    private TextField workplaceField;
    private ComboBox<String> countryComboBox;
    private TextField licenseKeyField;
    
    private Button registerButton;
    private Button activateButton;
    
    private RegistrationManager registrationManager;
    
    /**
     * Create a new registration dialog
     */
    public RegistrationDialog() {
        super();
        
        // Initialize registration manager
        registrationManager = RegistrationManager.getInstance();
        
        // Set up the dialog
        setTitle("MetaPilot Registration");
        setHeaderText("Registration is required to use MetaPilot");
        
        // Initialize components
        initComponents();
        
        // Set the result converter
        setResultConverter(buttonType -> {
            if (buttonType == ButtonType.OK) {
                return new UserRegistration(
                    nameField.getText(),
                    emailField.getText(),
                    workplaceField.getText(),
                    countryComboBox.getValue(),
                    licenseKeyField.getText()
                );
            }
            return null;
        });
        
        // Set dialog to be non-closable by the X button
        Stage stage = (Stage) getDialogPane().getScene().getWindow();
        stage.setOnCloseRequest(event -> {
            if (!isActivated()) {
                event.consume();
                Alert alert = new Alert(Alert.AlertType.WARNING);
                alert.setTitle("Registration Required");
                alert.setHeaderText("Registration is required to use MetaPilot");
                alert.setContentText("You must complete the registration process to use MetaPilot.");
                alert.showAndWait();
            }
        });
    }
    
    /**
     * Initialize the dialog components
     */
    private void initComponents() {
        // Create fields
        nameField = new TextField();
        nameField.setPromptText("Full Name");
        
        emailField = new TextField();
        emailField.setPromptText("Email Address");
        
        workplaceField = new TextField();
        workplaceField.setPromptText("Institution/Company");
        
        countryComboBox = new ComboBox<>();
        populateCountries();
        
        licenseKeyField = new TextField();
        licenseKeyField.setPromptText("License Key");
        licenseKeyField.setDisable(true);
        
        // Create register button
        registerButton = new Button("Register");
        registerButton.setDefaultButton(true);
        registerButton.setOnAction(e -> register());
        
        // Create activate button
        activateButton = new Button("Activate");
        activateButton.setDisable(true);
        activateButton.setOnAction(e -> activate());
        
        // Create the layout
        GridPane grid = new GridPane();
        grid.setHgap(10);
        grid.setVgap(10);
        grid.setPadding(new Insets(20, 150, 10, 10));
        
        // Add components to the grid
        grid.add(new Label("Name:"), 0, 0);
        grid.add(nameField, 1, 0);
        
        grid.add(new Label("Email:"), 0, 1);
        grid.add(emailField, 1, 1);
        
        grid.add(new Label("Workplace:"), 0, 2);
        grid.add(workplaceField, 1, 2);
        
        grid.add(new Label("Country:"), 0, 3);
        grid.add(countryComboBox, 1, 3);
        
        VBox buttonBox = new VBox(10);
        buttonBox.setAlignment(Pos.CENTER);
        buttonBox.getChildren().addAll(registerButton, activateButton);
        
        grid.add(buttonBox, 1, 4);
        
        grid.add(new Label("License Key:"), 0, 5);
        grid.add(licenseKeyField, 1, 5);
        
        // Create a section for terms and conditions
        CheckBox termsCheckBox = new CheckBox("I agree to the terms and conditions");
        Hyperlink viewTermsLink = new Hyperlink("View Terms and Conditions");
        viewTermsLink.setOnAction(e -> openTermsAndConditions());
        
        VBox termsBox = new VBox(5);
        termsBox.getChildren().addAll(termsCheckBox, viewTermsLink);
        grid.add(termsBox, 1, 6);
        
        // Create a registration instructions label
        Label instructionsLabel = new Label(
            "To use MetaPilot, please register by filling out the form above.\n" +
            "A license key will be generated when you click Register.\n" +
            "Click Activate to complete the registration process."
        );
        instructionsLabel.setWrapText(true);
        grid.add(instructionsLabel, 0, 7, 2, 1);
        
        // Enable/disable register button based on terms acceptance
        registerButton.disableProperty().bind(termsCheckBox.selectedProperty().not());
        
        // Add the grid to the dialog
        getDialogPane().setContent(grid);
        
        // Add buttons to the dialog
        getDialogPane().getButtonTypes().addAll(ButtonType.OK, ButtonType.CANCEL);
        
        // Hide the OK button initially - it will be enabled after activation
        Node okButton = getDialogPane().lookupButton(ButtonType.OK);
        okButton.setVisible(false);
        okButton.setManaged(false);
        
        // Set initial focus to the name field
        Platform.runLater(() -> nameField.requestFocus());
    }
    
    /**
     * Populate the country dropdown with a list of countries
     */
    private void populateCountries() {
        List<String> countries = Arrays.asList(Locale.getISOCountries());
        
        // Convert country codes to full names and sort alphabetically
        countries.stream()
            .map(code -> new Locale("", code).getDisplayCountry())
            .sorted()
            .forEach(countryComboBox.getItems()::add);
            
        // Set default selection to user's locale if available
        String userCountry = Locale.getDefault().getDisplayCountry();
        if (countryComboBox.getItems().contains(userCountry)) {
            countryComboBox.setValue(userCountry);
        } else {
            countryComboBox.getSelectionModel().selectFirst();
        }
    }
    
    /**
     * Open terms and conditions in the default browser
     */
    private void openTermsAndConditions() {
        try {
            // In a real application, this would point to your actual terms page
            // For now, we'll use a placeholder URL
            URI uri = new URI("https://example.com/metapilot/terms");
            Desktop.getDesktop().browse(uri);
        } catch (Exception e) {
            Alert alert = new Alert(Alert.AlertType.ERROR);
            alert.setTitle("Cannot Open Browser");
            alert.setHeaderText("Failed to open terms and conditions");
            alert.setContentText("Please visit our website manually to view the terms and conditions.");
            alert.showAndWait();
        }
    }
    
    /**
     * Handle the registration process
     */
    private void register() {
        // Validate input
        if (!validateInput()) {
            return;
        }
        
        // Create registration object
        UserRegistration registration = new UserRegistration(
            nameField.getText(),
            emailField.getText(),
            workplaceField.getText(),
            countryComboBox.getValue()
        );
        
        // Save to database
        boolean success = registrationManager.register(registration);
        
        if (success) {
            // Show license key
            licenseKeyField.setText(registration.getLicenseKey());
            licenseKeyField.setDisable(false);
            
            // Enable activation button
            activateButton.setDisable(false);
            
            // Show success message
            Alert alert = new Alert(Alert.AlertType.INFORMATION);
            alert.setTitle("Registration Successful");
            alert.setHeaderText("Thank you for registering MetaPilot");
            alert.setContentText("Your license key has been generated. Please use this key to activate the software.");
            alert.showAndWait();
        } else {
            // Show error message
            Alert alert = new Alert(Alert.AlertType.ERROR);
            alert.setTitle("Registration Failed");
            alert.setHeaderText("Failed to register");
            alert.setContentText("The email address may already be registered. Please try again or contact support.");
            alert.showAndWait();
        }
    }
    
    /**
     * Handle the activation process
     */
    private void activate() {
        String email = emailField.getText();
        String licenseKey = licenseKeyField.getText();
        
        // Attempt to activate
        boolean success = registrationManager.activateRegistration(email, licenseKey);
        
        if (success) {
            // Show success message
            Alert alert = new Alert(Alert.AlertType.INFORMATION);
            alert.setTitle("Activation Successful");
            alert.setHeaderText("MetaPilot has been activated");
            alert.setContentText("Thank you for registering. You can now use MetaPilot.");
            alert.showAndWait();
            
            // Enable OK button
            Node okButton = getDialogPane().lookupButton(ButtonType.OK);
            okButton.setVisible(true);
            okButton.setManaged(true);
            
            // Disable registration controls
            nameField.setDisable(true);
            emailField.setDisable(true);
            workplaceField.setDisable(true);
            countryComboBox.setDisable(true);
            registerButton.setDisable(true);
            activateButton.setDisable(true);
        } else {
            // Show error message
            Alert alert = new Alert(Alert.AlertType.ERROR);
            alert.setTitle("Activation Failed");
            alert.setHeaderText("Failed to activate MetaPilot");
            alert.setContentText("The license key may be invalid. Please check and try again.");
            alert.showAndWait();
        }
    }
    
    /**
     * Validate the registration input fields
     */
    private boolean validateInput() {
        StringBuilder errors = new StringBuilder();
        
        if (nameField.getText().trim().isEmpty()) {
            errors.append("- Name is required\n");
        }
        
        if (emailField.getText().trim().isEmpty()) {
            errors.append("- Email is required\n");
        } else if (!isValidEmail(emailField.getText().trim())) {
            errors.append("- Email format is invalid\n");
        }
        
        if (workplaceField.getText().trim().isEmpty()) {
            errors.append("- Workplace is required\n");
        }
        
        if (countryComboBox.getValue() == null) {
            errors.append("- Country is required\n");
        }
        
        if (errors.length() > 0) {
            Alert alert = new Alert(Alert.AlertType.ERROR);
            alert.setTitle("Validation Error");
            alert.setHeaderText("Please correct the following errors:");
            alert.setContentText(errors.toString());
            alert.showAndWait();
            return false;
        }
        
        return true;
    }
    
    /**
     * Check if the email format is valid
     */
    private boolean isValidEmail(String email) {
        // Simple email validation
        String regex = "^[A-Za-z0-9+_.-]+@(.+)$";
        return email.matches(regex);
    }
    
    /**
     * Check if registration is already activated
     */
    public boolean isActivated() {
        return licenseKeyField.isDisabled() && 
               nameField.isDisabled() && 
               emailField.isDisabled() && 
               workplaceField.isDisabled() && 
               countryComboBox.isDisabled();
    }
    
    /**
     * Show an existing registration
     */
    public void showExistingRegistration(String email, String licenseKey) {
        emailField.setText(email);
        licenseKeyField.setText(licenseKey);
        
        // Enable the activate button for re-activation
        activateButton.setDisable(false);
        licenseKeyField.setDisable(false);
    }
} 