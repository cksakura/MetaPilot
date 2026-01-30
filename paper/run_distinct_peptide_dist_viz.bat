@echo off
REM Run the distinct peptide distribution violin plot visualization.
REM Usage: run_distinct_peptide_dist_viz.bat
REM Set Python environment path
set PYTHON_PATH=C:\Users\daz25hak\MetaProAI\Resources\python-3.13.7-embed-amd64

set INPUT=distinct_peptide_distribution.tsv
set OUTPUT=distinct_peptide_distribution.html

"%PYTHON_PATH%\python.exe" "%~dp0distinct_peptide_dist_viz.py" "%~dp0%INPUT%" "%~dp0%OUTPUT%"

if %ERRORLEVEL% neq 0 (
    echo Visualization failed. Ensure Python, plotly, and pandas are installed.
    exit /b %ERRORLEVEL%
)

echo Violin plot written to %OUTPUT%

