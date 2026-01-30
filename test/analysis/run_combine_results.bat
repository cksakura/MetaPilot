@echo off
REM Combine result files from Figure1D, Figure1E, and Figure1F folders

set JAVA_HOME=C:\Users\daz25hak\AppData\Local\Programs\Eclipse Adoptium\jdk-21.0.8.9-hotspot
set PATH=%JAVA_HOME%\bin;%PATH%

echo.
echo ========================================
echo Combine Results
echo ========================================
echo.

cd /d "%~dp0"

echo Compiling CombineResults.java...
javac CombineResults.java

if not exist "CombineResults.class" (
    echo COMPILATION FAILED!
    pause
    exit /b 1
)

echo Compilation successful!
echo.
echo Combining files...
echo ========================================
echo.

java CombineResults

if %ERRORLEVEL% NEQ 0 (
    echo Combination failed!
    pause
    exit /b %ERRORLEVEL%
)

echo.
echo ========================================
echo Combination complete!
echo ========================================
echo.
pause

