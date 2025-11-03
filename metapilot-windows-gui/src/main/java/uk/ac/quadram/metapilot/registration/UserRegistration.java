package uk.ac.quadram.metapilot.registration;

/**
 * Data class for user registration information
 */
public class UserRegistration {
    private String name;
    private String email;
    private String workplace;
    private String country;
    private String licenseKey;
    
    public UserRegistration(String name, String email, String workplace, String country) {
        this.name = name;
        this.email = email;
        this.workplace = workplace;
        this.country = country;
        this.licenseKey = RegistrationManager.generateLicenseKey(email);
    }
    
    public UserRegistration(String name, String email, String workplace, String country, String licenseKey) {
        this.name = name;
        this.email = email;
        this.workplace = workplace;
        this.country = country;
        this.licenseKey = licenseKey;
    }
    
    public String getName() {
        return name;
    }
    
    public void setName(String name) {
        this.name = name;
    }
    
    public String getEmail() {
        return email;
    }
    
    public void setEmail(String email) {
        this.email = email;
    }
    
    public String getWorkplace() {
        return workplace;
    }
    
    public void setWorkplace(String workplace) {
        this.workplace = workplace;
    }
    
    public String getCountry() {
        return country;
    }
    
    public void setCountry(String country) {
        this.country = country;
    }
    
    public String getLicenseKey() {
        return licenseKey;
    }
    
    public void setLicenseKey(String licenseKey) {
        this.licenseKey = licenseKey;
    }
} 