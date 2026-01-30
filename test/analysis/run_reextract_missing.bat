@echo off
REM Re-extract .raw files for missing projects and append to existing list

set JAVA_HOME=C:\Users\daz25hak\AppData\Local\Programs\Eclipse Adoptium\jdk-21.0.8.9-hotspot
set PATH=%JAVA_HOME%\bin;%PATH%

echo.
echo ========================================
echo ReExtractMissingProjects - Re-extract Missing Projects
echo ========================================
echo.

cd /d "%~dp0"

echo Compiling ReExtractMissingProjects.java...
javac ReExtractMissingProjects.java

if not exist "ReExtractMissingProjects.class" (
    echo COMPILATION FAILED!
    pause
    exit /b 1
)

echo Compilation successful!
echo.
echo Running ReExtractMissingProjects...
echo ========================================
echo.

java ReExtractMissingProjects

if %ERRORLEVEL% NEQ 0 (
    echo.
    echo EXECUTION FAILED!
    pause
    exit /b %ERRORLEVEL%
)

echo.
echo Done! Check raw_files_list_complete.txt for the complete list.
pause

