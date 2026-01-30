@echo off
REM Run the protein and peptide distribution violin plot visualization.
REM Usage: run_protein_peptide_dist_viz.bat
REM Set Python environment path
set PYTHON_PATH=C:\Users\daz25hak\MetaProAI\Resources\python-3.13.7-embed-amd64

set INPUT=protein_peptide_distribution.tsv
set OUTPUT=protein_peptide_distribution.html

"%PYTHON_PATH%\python.exe" "%~dp0protein_peptide_dist_viz.py" "%~dp0%INPUT%" "%~dp0%OUTPUT%"

if %ERRORLEVEL% neq 0 (
    echo Visualization failed. Ensure Python, plotly, and pandas are installed.
    exit /b %ERRORLEVEL%
)

echo Violin plots written to %OUTPUT%

