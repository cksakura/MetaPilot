@echo off
REM Compile and run CountRarContents.java
REM Counts top-level entries in RAR files without downloading the entire archive

set JAVA_HOME=C:\Users\daz25hak\AppData\Local\Programs\Eclipse Adoptium\jdk-21.0.8.9-hotspot
set PATH=%JAVA_HOME%\bin;%PATH%

echo.
echo ========================================
echo CountRarContents - Count RAR File Contents
echo ========================================
echo.

cd /d "%~dp0"

echo Compiling CountRarContents.java...
javac CountRarContents.java

if not exist "CountRarContents.class" (
    echo COMPILATION FAILED!
    pause
    exit /b 1
)

echo Compilation successful!
echo.
echo Usage: java CountRarContents ^<rar_file_url^>
echo Example: java CountRarContents "https://ftp.pride.ebi.ac.uk/pride/data/archive/2016/08/PXD003791/file.rar"
echo.

if "%~1"=="" (
    echo No URL provided. Please provide a RAR file URL as argument.
    echo.
    echo Usage: run_count_rar.bat "https://url/to/file.rar"
    pause
    exit /b 1
)

echo Running CountRarContents...
echo ========================================
echo.

java CountRarContents "%~1"

pause

