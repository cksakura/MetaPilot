package uk.ac.quadram.metapilot.registration;

import java.io.File;
import java.sql.Connection;
import java.sql.DriverManager;
import java.sql.PreparedStatement;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.sql.Statement;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * Manages user registration and storage in a SQLite database
 */
public class RegistrationManager {
    private static final Logger LOGGER = Logger.getLogger(RegistrationManager.class.getName());
    private static final String DB_FOLDER = "config";
    private static final String DB_NAME = "registration.db";
    private static RegistrationManager instance;
    
    private String dbPath;
    private Connection connection;
    
    /**
     * Get singleton instance of the registration manager
     */
    public static synchronized RegistrationManager getInstance() {
        if (instance == null) {
            instance = new RegistrationManager();
        }
        return instance;
    }
    
    /**
     * Private constructor for singleton pattern
     */
    private RegistrationManager() {
        // Initialize the database
        initializeDatabase();
    }
    
    /**
     * Initialize the SQLite database
     */
    private void initializeDatabase() {
        try {
            // Ensure config directory exists
            File configDir = new File(DB_FOLDER);
            if (!configDir.exists()) {
                configDir.mkdirs();
            }
            
            // Set database path
            dbPath = DB_FOLDER + File.separator + DB_NAME;
            
            // Load JDBC driver
            Class.forName("org.sqlite.JDBC");
            
            // Get connection
            connection = DriverManager.getConnection("jdbc:sqlite:" + dbPath);
            
            // Create tables if they don't exist
            createTables();
            
            LOGGER.info("Registration database initialized at: " + dbPath);
        } catch (ClassNotFoundException | SQLException e) {
            LOGGER.log(Level.SEVERE, "Failed to initialize registration database", e);
        }
    }
    
    /**
     * Create required database tables
     */
    private void createTables() throws SQLException {
        try (Statement statement = connection.createStatement()) {
            // Create registration table
            statement.execute(
                "CREATE TABLE IF NOT EXISTS registration (" +
                "id INTEGER PRIMARY KEY AUTOINCREMENT," +
                "name TEXT NOT NULL," +
                "email TEXT NOT NULL UNIQUE," +
                "workplace TEXT NOT NULL," +
                "country TEXT NOT NULL," +
                "registration_date TIMESTAMP DEFAULT CURRENT_TIMESTAMP," +
                "license_key TEXT NOT NULL UNIQUE" +
                ");"
            );
            
            // Create activation table
            statement.execute(
                "CREATE TABLE IF NOT EXISTS activation (" +
                "id INTEGER PRIMARY KEY AUTOINCREMENT," +
                "email TEXT NOT NULL UNIQUE," +
                "activation_date TIMESTAMP DEFAULT CURRENT_TIMESTAMP," +
                "FOREIGN KEY (email) REFERENCES registration(email)" +
                ");"
            );
        }
    }
    
    /**
     * Save a new registration to the database
     * @return true if registration was successful
     */
    public boolean register(UserRegistration registration) {
        String sql = "INSERT INTO registration (name, email, workplace, country, license_key) VALUES (?, ?, ?, ?, ?)";
        
        try (PreparedStatement stmt = connection.prepareStatement(sql)) {
            stmt.setString(1, registration.getName());
            stmt.setString(2, registration.getEmail());
            stmt.setString(3, registration.getWorkplace());
            stmt.setString(4, registration.getCountry());
            stmt.setString(5, registration.getLicenseKey());
            
            int rowsAffected = stmt.executeUpdate();
            return rowsAffected > 0;
        } catch (SQLException e) {
            LOGGER.log(Level.SEVERE, "Error registering user", e);
            return false;
        }
    }
    
    /**
     * Activate a registration with the given email and license key
     * @return true if activation was successful
     */
    public boolean activateRegistration(String email, String licenseKey) {
        // First verify the license key
        if (!verifyLicenseKey(email, licenseKey)) {
            return false;
        }
        
        // Then record the activation
        String sql = "INSERT OR REPLACE INTO activation (email) VALUES (?)";
        
        try (PreparedStatement stmt = connection.prepareStatement(sql)) {
            stmt.setString(1, email);
            int rowsAffected = stmt.executeUpdate();
            return rowsAffected > 0;
        } catch (SQLException e) {
            LOGGER.log(Level.SEVERE, "Error activating registration", e);
            return false;
        }
    }
    
    /**
     * Verify if a license key is valid for the given email
     */
    public boolean verifyLicenseKey(String email, String licenseKey) {
        String sql = "SELECT COUNT(*) FROM registration WHERE email = ? AND license_key = ?";
        
        try (PreparedStatement stmt = connection.prepareStatement(sql)) {
            stmt.setString(1, email);
            stmt.setString(2, licenseKey);
            
            try (ResultSet rs = stmt.executeQuery()) {
                if (rs.next()) {
                    return rs.getInt(1) > 0;
                }
            }
        } catch (SQLException e) {
            LOGGER.log(Level.SEVERE, "Error verifying license key", e);
        }
        
        return false;
    }
    
    /**
     * Check if a user is registered and activated
     */
    public boolean isRegisteredAndActivated(String email) {
        String sql = "SELECT COUNT(*) FROM activation WHERE email = ?";
        
        try (PreparedStatement stmt = connection.prepareStatement(sql)) {
            stmt.setString(1, email);
            
            try (ResultSet rs = stmt.executeQuery()) {
                if (rs.next()) {
                    return rs.getInt(1) > 0;
                }
            }
        } catch (SQLException e) {
            LOGGER.log(Level.SEVERE, "Error checking registration status", e);
        }
        
        return false;
    }
    
    /**
     * Generate a unique license key for a user
     * This is a simple implementation that should be made more secure in production
     */
    public static String generateLicenseKey(String email) {
        // Simple hashing of email with a timestamp to create a license key
        // In a real application, you'd want a more secure algorithm
        String timestamp = String.valueOf(System.currentTimeMillis());
        String combined = email + timestamp;
        
        // Create a simple hash
        int hash = combined.hashCode();
        String licenseKey = String.format("MP-%08X-%04X-%04X", 
                hash & 0xFFFFFFFF, 
                (hash >> 16) & 0xFFFF,
                (int)(Math.random() * 0xFFFF));
        
        return licenseKey;
    }
    
    /**
     * Close the database connection
     */
    public void close() {
        try {
            if (connection != null && !connection.isClosed()) {
                connection.close();
            }
        } catch (SQLException e) {
            LOGGER.log(Level.SEVERE, "Error closing database connection", e);
        }
    }
} 