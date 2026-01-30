@echo off
REM Extract project information from pride_datasets.json to JSON and CSV formats

set JAVA_HOME=C:\Users\daz25hak\AppData\Local\Programs\Eclipse Adoptium\jdk-21.0.8.9-hotspot
set PATH=%JAVA_HOME%\bin;%PATH%

echo.
echo ========================================
echo Extract Project Information
echo ========================================
echo.

cd /d "%~dp0"

echo Compiling ExtractProjectInfo.java...
javac ExtractProjectInfo.java

if not exist "ExtractProjectInfo.class" (
    echo COMPILATION FAILED!
    pause
    exit /b 1
)

echo Compilation successful!
echo.
echo Extracting to JSON format...
echo ========================================
echo.

java ExtractProjectInfo

if %ERRORLEVEL% NEQ 0 (
    echo JSON extraction failed!
    pause
    exit /b %ERRORLEVEL%
)

echo.
echo Compiling ExtractProjectInfoCSV.java...
javac ExtractProjectInfoCSV.java

if not exist "ExtractProjectInfoCSV.class" (
    echo CSV compilation failed, but JSON extraction succeeded!
    pause
    exit /b 1
)

echo.
echo Extracting to CSV format...
echo ========================================
echo.

java ExtractProjectInfoCSV

if %ERRORLEVEL% NEQ 0 (
    echo CSV extraction failed, but JSON extraction succeeded!
    pause
    exit /b %ERRORLEVEL%
)

echo.
echo ========================================
echo Extraction complete!
echo ========================================
echo.
echo Output files:
echo   - selected_projects_info.json (JSON format)
echo   - selected_projects_info.csv (CSV format for Excel/Python)
echo.
pause

