@echo off
echo Starting MetaPilot Application...

:: Set Java options
set JAVA_OPTS=-Xmx2g

:: Run the application with the fat jar that contains all dependencies
java %JAVA_OPTS% --module-path "C:\Program Files\Java\javafx-sdk-21.0.7\lib" --add-modules javafx.controls,javafx.fxml,javafx.graphics,javafx.web,javafx.swing,javafx.base -jar target/metapilot-1.0-SNAPSHOT.jar

pause 