@echo off
REM Compile and run DatasetMarkerAnalysis for experimental dataset analysis
REM This analyzes peptide-intensity files to calculate marker-4 protein abundance ratios

set JAVA_HOME=C:\Users\daz25hak\AppData\Local\Programs\Eclipse Adoptium\jdk-21.0.8.9-hotspot
set PATH=%JAVA_HOME%\bin;%PATH%

echo Compiling DatasetMarkerAnalysis and MarkerGenomeSelector...
javac -d . ..\..\metapilot-core\src\main\java\uk\ac\quadram\metapilot\task\mag\MarkerGenomeSelector.java DatasetMarkerAnalysis.java

if %ERRORLEVEL% neq 0 (
    echo Compilation failed!
    exit /b %ERRORLEVEL%
)

echo Running dataset analysis...
java analysis.DatasetMarkerAnalysis

if %ERRORLEVEL% neq 0 (
    echo Analysis failed!
    exit /b %ERRORLEVEL%
)

echo Analysis complete! Check output files

