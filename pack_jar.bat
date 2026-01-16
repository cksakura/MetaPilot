@echo off
REM ===================================================================
REM MetaPilot JAR Packing Script
REM 
REM This script builds and packages the MetaPilot Windows GUI JAR file.
REM
REM Usage:
REM   pack_jar.bat        Build JAR (with clean)
REM ===================================================================

setlocal enabledelayedexpansion

echo ================================================
echo   MetaPilot JAR Packing
echo ================================================
echo.

REM Check Java
echo Checking Java...
java -version >nul 2>&1
if errorlevel 1 (
    echo ERROR: Java is not installed or not in PATH
    pause
    exit /b 1
)
echo [OK] Java found

REM Check Maven
echo Checking Maven...
mvn -version >nul 2>&1
if errorlevel 1 (
    echo ERROR: Maven is not installed or not in PATH
    pause
    exit /b 1
)
echo [OK] Maven found
echo.

REM Build
echo Building MetaPilot Windows GUI JAR...
echo This may take a few minutes...
echo.

REM Set MAVEN_OPTS for Maven JVM memory
set MAVEN_OPTS=-Xmx2g

REM Build with clean and package, skip tests, build dependencies
call mvn clean package -am -DskipTests -pl metapilot-windows-gui

if errorlevel 1 (
    echo.
    echo ERROR: Build failed
    pause
    exit /b 1
)

echo.
echo ================================================
echo   Build Complete!
echo ================================================
echo.

REM Check for the JAR file
set JAR_FILE=metapilot-windows-gui\target\metapilot-windows-gui-1.0.0-beta.5-shaded.jar
if exist "%JAR_FILE%" (
    echo [OK] JAR file created successfully!
    echo.
    echo Location: %JAR_FILE%
    echo.
    echo To run the application:
    echo    java -jar "%JAR_FILE%"
    echo.
) else (
    echo WARNING: Expected JAR file not found at: %JAR_FILE%
    echo Please check the build output above for errors.
    echo.
)

pause
endlocal
