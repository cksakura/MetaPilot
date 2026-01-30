@echo off
REM Compile and run ProjectStatistics.java
REM Analyzes raw_files_list.txt to generate per-year statistics

set JAVA_HOME=C:\Users\daz25hak\AppData\Local\Programs\Eclipse Adoptium\jdk-21.0.8.9-hotspot
set PATH=%JAVA_HOME%\bin;%PATH%

echo.
echo ========================================
echo ProjectStatistics - Yearly Statistics
echo ========================================
echo.

REM Change to script directory
cd /d "%~dp0"

echo Compiling ProjectStatistics.java...
javac ProjectStatistics.java

if not exist "ProjectStatistics.class" (
    echo.
    echo COMPILATION FAILED!
    echo.
    pause
    exit /b 1
)

echo Compilation successful!
echo.
echo Running ProjectStatistics...
echo ========================================
echo.

java ProjectStatistics

if %ERRORLEVEL% NEQ 0 (
    echo.
    echo ========================================
    echo EXECUTION FAILED with exit code %ERRORLEVEL%
    echo ========================================
    pause
    exit /b %ERRORLEVEL%
)

echo.
echo ========================================
echo Done! Check project_statistics.txt for results.
echo ========================================
pause

