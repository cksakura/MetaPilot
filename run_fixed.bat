@echo off
echo Starting MetaPilot Application with Debug Options...

set JAVA_OPTS=-Xmx2g -XX:+ShowCodeDetailsInExceptionMessages -Dlog4j.configurationFile=file:src/main/resources/log4j2.xml -Dlog4j2.debug=true

java %JAVA_OPTS% --module-path "C:\Program Files\Java\javafx-sdk-21.0.7\lib" --add-modules javafx.controls,javafx.fxml,javafx.graphics,javafx.web,javafx.swing,javafx.base -cp "target/classes;lib/*;lib/miglayout-swing-5.3.jar;lib/miglayout-core-5.3.jar" uk.ac.quadram.metapilot.MetaPilotApp

pause 