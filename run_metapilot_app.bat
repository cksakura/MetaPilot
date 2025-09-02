
@echo off
setlocal
echo Starting MetaPilot Application...

:: Set Java and Maven paths
set "JAVA_HOME=C:\Users\daz25hak\AppData\Local\Programs\Eclipse Adoptium\jdk-21.0.8.9-hotspot"
set "MAVEN_HOME=C:\Program Files\Apache\apache-maven-3.9.11-bin\apache-maven-3.9.11"
set "PATH=%JAVA_HOME%\bin;%MAVEN_HOME%\bin;%PATH%"

:: Set Java options (override by setting JAVA_OPTS before calling)
if not defined JAVA_OPTS set "JAVA_OPTS=-Xmx2g --add-opens java.xml/javax.xml.parsers=ALL-UNNAMED --add-opens java.xml/org.w3c.dom=ALL-UNNAMED --add-opens java.desktop/java.awt=ALL-UNNAMED --add-opens java.base/java.lang=ALL-UNNAMED"

:: JavaFX path (override by setting JAVAFX_DIR before calling)
if not defined JAVAFX_DIR set "JAVAFX_DIR=C:\Program Files\Java\openjfx-21.0.8_windows-x64_bin-sdk\javafx-sdk-21.0.8\lib"

:: Shaded JAR produced by Maven Shade Plugin
set "JAR=target\metapilot-1.0-SNAPSHOT.jar"

:: Check if Java is available
java -version >nul 2>&1
if errorlevel 1 (
  echo ERROR: Java is not installed or not in PATH.
  echo Please install Java 17 or later from: https://adoptium.net/
  echo After installation, restart your terminal and try again.
  if "%NO_PAUSE%"=="" pause
  exit /b 1
)

:: Check if JAR exists
if not exist "%JAR%" (
  echo JAR not found at "%JAR%". Attempting to build with Maven...
  
  :: Check if Maven is available using full path first
  "%MAVEN_HOME%\bin\mvn.cmd" -version >nul 2>&1
  if errorlevel 1 (
    echo ERROR: Maven is not available at "%MAVEN_HOME%\bin\mvn.cmd"
    echo.
    echo Please verify that Maven is installed at:
    echo %MAVEN_HOME%
    echo.
    echo Or update the MAVEN_HOME path in this batch file.
    echo.
    echo Alternatively, if you have a pre-built JAR file, place it at:
    echo %JAR%
    if "%NO_PAUSE%"=="" pause
    exit /b 1
  )
  
  :: Build with Maven using full path
  echo Building project with Maven...
  "%MAVEN_HOME%\bin\mvn.cmd" -q -DskipTests package
  if errorlevel 1 (
    echo ERROR: Maven build failed.
    if "%NO_PAUSE%"=="" pause
    exit /b 1
  )
  
  :: Verify JAR exists after build
  if not exist "%JAR%" (
    echo ERROR: Build completed but JAR not found at "%JAR%".
    echo Check the Maven output for errors.
    if "%NO_PAUSE%"=="" pause
    exit /b 1
  )
)

:: Run the main class directly
echo Running MetaPilot main class: uk.ac.quadram.metapilot.MetaPilotApp

:: Try running with JavaFX first
if exist "%JAVAFX_DIR%" (
  echo Using JavaFX from: %JAVAFX_DIR%
  java %JAVA_OPTS% --module-path "%JAVAFX_DIR%" --add-modules javafx.controls,javafx.fxml,javafx.graphics,javafx.web,javafx.swing,javafx.base -cp "%JAR%" uk.ac.quadram.metapilot.MetaPilotApp
  if errorlevel 1 (
    echo.
    echo JavaFX version compatibility issue detected.
    echo Attempting to run without JavaFX modules...
    echo.
    java %JAVA_OPTS% -cp "%JAR%" uk.ac.quadram.metapilot.MetaPilotApp
  )
) else (
  echo WARNING: JavaFX directory not found at: %JAVAFX_DIR%
  echo You may need to download JavaFX SDK from: https://openjfx.io/
  echo Or set JAVAFX_DIR environment variable to the correct path.
  echo.
  echo Attempting to run without JavaFX modules...
  java %JAVA_OPTS% -cp "%JAR%" uk.ac.quadram.metapilot.MetaPilotApp
)

set "EXIT_CODE=%ERRORLEVEL%"
if not "%EXIT_CODE%"=="0" (
  echo Application exited with code %EXIT_CODE%.
)

if "%NO_PAUSE%"=="" pause
endlocal & exit /b %EXIT_CODE% 