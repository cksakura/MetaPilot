@echo off
REM Compile and run PeptideCategoryAnalyzer for paper analysis
REM This is standalone code, separate from main MetaPilot software

set JAVA_HOME=C:\Users\daz25hak\AppData\Local\Programs\Eclipse Adoptium\jdk-21.0.8.9-hotspot
set PATH=%JAVA_HOME%\bin;%PATH%

echo Compiling PeptideCategoryAnalyzer...
javac -d . PeptideCategoryAnalyzer.java

if %ERRORLEVEL% neq 0 (
    echo Compilation failed!
    exit /b %ERRORLEVEL%
)

echo Running analysis...
java analysis.PeptideCategoryAnalyzer

if %ERRORLEVEL% neq 0 (
    echo Analysis failed!
    exit /b %ERRORLEVEL%
)

echo Analysis complete! Check humangut_4category_peptides.csv

