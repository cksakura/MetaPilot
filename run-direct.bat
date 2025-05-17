@echo off
set JAVA_HOME="C:\Program Files\Java\jdk-17"
set JAVAFX_PATH="C:\Program Files\Java\javafx-sdk-21.0.7\lib"
set PATH=%JAVA_HOME%\bin;%PATH%

%JAVA_HOME%\bin\java --module-path %JAVAFX_PATH% --add-modules javafx.controls,javafx.fxml,javafx.web,javafx.media,javafx.swing ^
  --add-opens javafx.controls/com.sun.javafx.scene.control.behavior=ALL-UNNAMED ^
  --add-opens javafx.controls/com.sun.javafx.scene.control=ALL-UNNAMED ^
  --add-opens javafx.graphics/com.sun.javafx.scene=ALL-UNNAMED ^
  -jar target/MetaPilot-1.0-SNAPSHOT.jar 