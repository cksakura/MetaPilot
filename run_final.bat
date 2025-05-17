@echo off
echo Starting MetaPilot Application...

:: Set Java options
set JAVA_OPTS=-Xmx2g -Dlog4j.configurationFile=file:src/main/resources/log4j2.xml

:: Ensure logs directory exists
if not exist logs mkdir logs

:: Run the application
java %JAVA_OPTS% --module-path "C:\Program Files\Java\javafx-sdk-21.0.7\lib" --add-modules javafx.controls,javafx.fxml,javafx.graphics,javafx.web,javafx.swing,javafx.base -cp "target/classes;lib/*" uk.ac.quadram.metapilot.MetaPilotApp

pause 