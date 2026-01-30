@echo off
REM Run the peptide sharing violin plot visualization (Plotly).
REM Usage: run_dist_pep_summary_violin.bat
REM Set Python environment path
set PYTHON_PATH=C:\Users\daz25hak\MetaProAI\Resources\python-3.13.7-embed-amd64

set INPUT=dist_pep_summary_humangut.tsv
set OUTPUT=dist_pep_summary_humangut_violin.html

"%PYTHON_PATH%\python.exe" "%~dp0dist_pep_summary_violin.py" "%~dp0%INPUT%" "%~dp0%OUTPUT%"

if %ERRORLEVEL% neq 0 (
    echo Visualization failed. Ensure Python, plotly, and pandas are installed.
    exit /b %ERRORLEVEL%
)

echo Violin plot written to %OUTPUT%

