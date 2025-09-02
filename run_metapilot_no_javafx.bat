@echo off
setlocal
echo Starting MetaPilot Application (without JavaFX)...

:: Set Java and Maven paths
set "JAVA_HOME=C:\Users\daz25hak\AppData\Local\Programs\Eclipse Adoptium\jdk-21.0.8.9-hotspot"
set "MAVEN_HOME=C:\Program Files\Apache\apache-maven-3.9.11-bin\apache-maven-3.9.11"
set "PATH=%JAVA_HOME%\bin;%MAVEN_HOME%\bin;%PATH%"

:: Set Java options (override by setting JAVA_OPTS before calling)
if not defined JAVA_OPTS set "JAVA_OPTS=-Xmx2g"

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

:: Run the main class without JavaFX
echo Running MetaPilot main class without JavaFX: uk.ac.quadram.metapilot.MetaPilotApp
echo Note: This may not work if the application requires JavaFX GUI components.
echo.
java %JAVA_OPTS% -cp "%JAR%" uk.ac.quadram.metapilot.MetaPilotApp

set "EXIT_CODE=%ERRORLEVEL%"
if not "%EXIT_CODE%"=="0" (
  echo Application exited with code %EXIT_CODE%.
  echo.
  echo If you see JavaFX-related errors, you may need to:
  echo 1. Download JavaFX 21.0.7 (compatible with Java 21) from: https://openjfx.io/
  echo 2. Or upgrade to Java 22 to use your current JavaFX 24.0.2
)

if "%NO_PAUSE%"=="" pause
endlocal & exit /b %EXIT_CODE%
