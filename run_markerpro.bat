@echo off
REM Batch file to run test_markerpro.py

REM Set Python environment path
set PYTHON_PATH=C:\Users\daz25hak\MetaProAI\Resources\python-3.13.7-embed-amd64

REM Prompt user for the TSV file path
echo ========================================
echo  MarkerPro Analysis Script
echo ========================================
echo.
set /p TSV_PATH="Enter the full path to functions.tsv file: "

REM Check if path was provided
if "%TSV_PATH%"=="" (
    echo Error: No file path provided.
    pause
    exit /b 1
)

REM Change to the script directory
cd /d "%~dp0src\main\resources\python_scripts"

REM Run the Python script with the provided path
echo.
echo Running analysis on: %TSV_PATH%
echo.
"%PYTHON_PATH%\python.exe" test_markerpro.py "%TSV_PATH%"

REM Pause to see output
echo.
echo Script execution completed.
pause

