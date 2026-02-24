@echo off
REM MetaPilot Build Script
REM Builds the MetaPilot Windows GUI module and generates the shaded JAR file.

setlocal enabledelayedexpansion

echo ========================================
echo Building MetaPilot
echo ========================================
echo.

REM Change to script directory (project root)
cd /d "%~dp0"

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

REM Set Maven JVM memory
set MAVEN_OPTS=-Xmx4g

REM Clean and build metapilot-windows-gui (and dependencies), skip tests
call mvn clean package -am -DskipTests -pl metapilot-windows-gui

REM Check if the build was successful
if %ERRORLEVEL% EQU 0 (
    echo.
    echo ========================================
    echo Build completed successfully!
    echo ========================================
    echo.
    echo The shaded JAR file has been generated at:
    echo metapilot-windows-gui\target\metapilot-windows-gui-*-shaded.jar
    echo.
    echo To run: java -jar metapilot-windows-gui\target\metapilot-windows-gui-*-shaded.jar
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
endlocal

