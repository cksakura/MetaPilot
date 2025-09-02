@echo off
setlocal
echo Starting MetaPilot Command Line Interface...

:: Set Java and Maven paths
set "JAVA_HOME=C:\Users\daz25hak\AppData\Local\Programs\Eclipse Adoptium\jdk-21.0.8.9-hotspot"
set "MAVEN_HOME=C:\Program Files\Apache\apache-maven-3.9.11-bin\apache-maven-3.9.11"
set "PATH=%JAVA_HOME%\bin;%MAVEN_HOME%\bin;%PATH%"

:: Set Java options (override by setting JAVA_OPTS before calling)
if not defined JAVA_OPTS set "JAVA_OPTS=-Xmx2g --add-opens java.xml/javax.xml.parsers=ALL-UNNAMED --add-opens java.xml/org.w3c.dom=ALL-UNNAMED --add-opens java.desktop/java.awt=ALL-UNNAMED --add-opens java.base/java.lang=ALL-UNNAMED"

:: JavaFX path for GUI option
set "JAVAFX_DIR=C:\Program Files\Java\openjfx-21.0.8_windows-x64_bin-sdk\javafx-sdk-21.0.8\lib"

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
  echo ERROR: JAR file not found at "%JAR%".
  echo Please run run_metapilot_app.bat first to build the project.
  if "%NO_PAUSE%"=="" pause
  exit /b 1
)

:: Run the command-line interface
if "%1"=="" (
  echo Usage: %0 [options]
  echo.
  echo Options:
  echo   --help     Show this help message
  echo   --gui      Launch the GUI (requires JavaFX)
  echo   --par file Run with parameter file
  echo   --sep file Run with separate parameter file
  echo.
  echo Example:
  echo   %0 --help
  echo   %0 --gui
  echo   %0 --par parameters.json
  echo.
  java %JAVA_OPTS% -cp "%JAR%" uk.ac.quadram.metapilot.task.v1.MetaLabCli --help
) else if "%1"=="--gui" (
  echo Launching GUI with JavaFX...
  if exist "%JAVAFX_DIR%" (
    java %JAVA_OPTS% --module-path "%JAVAFX_DIR%" --add-modules javafx.controls,javafx.fxml,javafx.graphics,javafx.web,javafx.swing,javafx.base -cp "%JAR%" uk.ac.quadram.metapilot.task.v1.MetaLabCli --gui
  ) else (
    echo WARNING: JavaFX not found, attempting to run without JavaFX...
    java %JAVA_OPTS% -cp "%JAR%" uk.ac.quadram.metapilot.task.v1.MetaLabCli --gui
  )
) else (
  java %JAVA_OPTS% -cp "%JAR%" uk.ac.quadram.metapilot.task.v1.MetaLabCli %*
)

set "EXIT_CODE=%ERRORLEVEL%"
if not "%EXIT_CODE%"=="0" (
  echo Application exited with code %EXIT_CODE%.
)

if "%NO_PAUSE%"=="" pause
endlocal & exit /b %EXIT_CODE%

