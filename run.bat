@echo off
set JAVAFX_PATH="C:\Program Files\Java\javafx-sdk-21.0.7\lib"
mvn clean javafx:run -Djavafx.modulepath="C:\Program Files\Java\javafx-sdk-21.0.7\lib" -Djavafx.module=javafx.controls