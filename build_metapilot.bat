@echo off
REM MetaPilot Build Script
REM This script builds the MetaPilot project and generates the shaded JAR file

echo ========================================
echo Building MetaPilot
echo ========================================
echo.

REM Check if Maven is available
where mvn >nul 2>nul
if %ERRORLEVEL% NEQ 0 (
    echo ERROR: Maven is not found in PATH
    echo Please install Maven or add it to your PATH environment variable
    pause
    exit /b 1
)

echo Maven found. Starting build...
echo.

REM Clean and build the project with Maven
mvn clean package

REM Check if the build was successful
if %ERRORLEVEL% EQU 0 (
    echo.
    echo ========================================
    echo Build completed successfully!
    echo ========================================
    echo.
    echo The shaded JAR file has been generated at:
    echo target\metapilot-1.0-SNAPSHOT-shaded.jar
    echo.
) else (
    echo.
    echo ========================================
    echo Build failed!
    echo ========================================
    echo Please check the error messages above.
    echo.
)

pause

