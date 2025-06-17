@echo off
setlocal

:: Set Java and Maven paths
set "JAVA_HOME=C:\Program Files\Java\jdk-21"
set "MAVEN_HOME=C:\Program Files\apache-maven-3.9.10"
set "PATH=%JAVA_HOME%\bin;%MAVEN_HOME%\bin;%PATH%"

:: Set JavaFX path (update this path to match your JavaFX SDK location)
set "JAVAFX_PATH=C:\Program Files\Java\javafx-sdk-21"

:: Display versions
echo Java version:
java -version
echo.
echo Maven version:
mvn -version
echo.

:: Run Maven build
echo Building project...
call mvn clean install

:: Check if build was successful
if %ERRORLEVEL% EQU 0 (
    echo.
    echo Build completed successfully!
    echo.
    echo Running application...
    java --module-path "%JAVAFX_PATH%\lib" --add-modules javafx.controls,javafx.fxml -cp target/metapilot-0.0.1-SNAPSHOT.jar uk.ac.quadram.metapilot.gui.MetaPilotApp
) else (
    echo.
    echo Build failed with error code %ERRORLEVEL%
)

pause 