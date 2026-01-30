@echo off
REM Compile and run AnalyzePrideDatasets.java
REM Analyzes PRIDE datasets JSON and generates summaries

set JAVA_HOME=C:\Users\daz25hak\AppData\Local\Programs\Eclipse Adoptium\jdk-21.0.8.9-hotspot
set PATH=%JAVA_HOME%\bin;%PATH%

echo.
echo ========================================
echo AnalyzePrideDatasets - PRIDE Dataset Analysis
echo ========================================
echo.

cd /d "%~dp0"

echo Compiling AnalyzePrideDatasets.java...
javac AnalyzePrideDatasets.java

if not exist "AnalyzePrideDatasets.class" (
    echo COMPILATION FAILED!
    pause
    exit /b 1
)

echo Compilation successful!
echo.
echo Running AnalyzePrideDatasets...
echo This may take a while for large JSON files...
echo ========================================
echo.

java AnalyzePrideDatasets

if %ERRORLEVEL% NEQ 0 (
    echo.
    echo EXECUTION FAILED!
    pause
    exit /b %ERRORLEVEL%
)

echo.
echo ========================================
echo Analysis complete! Check output files.
echo ========================================
echo.
echo To create word cloud visualization, run:
echo   python create_wordcloud.py
echo.
pause

