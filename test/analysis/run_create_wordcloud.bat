@echo off
REM Create word cloud visualization from word_cloud_data.csv

set PYTHON_PATH=C:\Users\daz25hak\AppData\Local\Programs\Python\Python312\python.exe

echo.
echo ========================================
echo Creating Word Cloud Visualization
echo ========================================
echo.

cd /d "%~dp0"

if not exist "word_cloud_data.csv" (
    echo Error: word_cloud_data.csv not found
    echo Please run AnalyzePrideDatasets.java first to generate the CSV file
    pause
    exit /b 1
)

echo Installing required packages (if needed)...
%PYTHON_PATH% -m pip install wordcloud matplotlib pandas --quiet

echo.
echo Creating word cloud...
%PYTHON_PATH% create_wordcloud.py

if %ERRORLEVEL% NEQ 0 (
    echo.
    echo FAILED to create word cloud!
    pause
    exit /b %ERRORLEVEL%
)

echo.
echo Word cloud created successfully!
echo Check wordcloud.png and wordcloud.html
pause

