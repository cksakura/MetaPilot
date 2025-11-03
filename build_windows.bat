@echo off
REM ===================================================================
REM MetaPilot Windows Build Script
REM 
REM This script builds the MetaPilot Windows GUI application.
REM
REM Usage:
REM   build_windows.bat        Build Windows GUI
REM   build_windows.bat clean  Clean before building
REM ===================================================================

setlocal enabledelayedexpansion

REM Configuration
set PROJECT_NAME=MetaPilot
set VERSION=1.0.0-beta.2
set BUILD_PROFILE=windows-gui

echo ================================================
echo   %PROJECT_NAME% Windows Build
echo   Version: %VERSION%
echo ================================================
echo.

REM Check Java
echo Checking requirements...
java -version >nul 2>&1
if errorlevel 1 (
    echo ERROR: Java is not installed or not in PATH
    exit /b 1
)
echo [OK] Java found

REM Check Maven
mvn -version >nul 2>&1
if errorlevel 1 (
    echo ERROR: Maven is not installed or not in PATH
    exit /b 1
)
echo [OK] Maven found
echo.

REM Build
set BUILD_TYPE=%1
if "%BUILD_TYPE%"=="clean" (
    echo Cleaning previous builds...
    call mvn clean -q
    echo [OK] Clean complete
    echo.
)

echo Building MetaPilot Windows GUI...
REM Set MAVEN_OPTS as environment variable for Maven JVM memory
set MAVEN_OPTS=-Xmx2g
REM Build with -am to also build dependencies (metapilot-core)
call mvn clean package -P%BUILD_PROFILE% -am -DskipTests -pl metapilot-windows-gui

if errorlevel 1 (
    echo.
    echo ERROR: Build failed
    exit /b 1
)

echo.
echo ================================================
echo   Build Complete!
echo ================================================
echo.
echo Build artifacts:
if exist "metapilot-windows-gui\target\metapilot-windows-gui-%VERSION%-shaded.jar" (
    echo [OK] Windows GUI JAR created
    echo    Location: metapilot-windows-gui\target\
    echo.
    echo To run:
    echo    java -jar metapilot-windows-gui\target\metapilot-windows-gui-%VERSION%-shaded.jar
    echo.
)

endlocal

