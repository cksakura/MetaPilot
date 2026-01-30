@echo off
REM Run the protein abundance statistical analysis.
REM Usage: run_protein_abundance_analysis.bat
REM Set Python environment path
set PYTHON_PATH=C:\Users\daz25hak\MetaProAI\Resources\python-3.13.7-embed-amd64

set INPUT=O480_60min_DIA_IGC_humanswiss_protein_abundances.tsv
set OUTPUT_PREFIX=O480_60min_DIA_IGC_humanswiss_protein_abundances_analysis

"%PYTHON_PATH%\python.exe" "%~dp0protein_abundance_analysis.py" "%~dp0%INPUT%" "%~dp0%OUTPUT_PREFIX%"

if %ERRORLEVEL% neq 0 (
    echo Analysis failed. Ensure Python, pandas, numpy, and plotly are installed.
    exit /b %ERRORLEVEL%
)

echo Analysis complete! Check the output files:
echo   - %OUTPUT_PREFIX%.stats.txt (statistical summary)
echo   - %OUTPUT_PREFIX%.html (interactive visualizations)

