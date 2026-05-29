@echo off
REM MetaPilot Build Script (CLI)
REM Builds the MetaPilot headless CLI module and generates the shaded JAR file.

setlocal enabledelayedexpansion

echo ========================================
echo Building MetaPilot CLI
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

REM Clean and build metapilot-cli (and dependencies), skip tests
call mvn clean package -am -DskipTests -pl metapilot-cli

REM Check if the build was successful
if %ERRORLEVEL% EQU 0 (
    echo.
    echo ========================================
    echo Build completed successfully!
    echo ========================================
    echo.
    REM Try to locate the newest shaded JAR
    set "SHADED_JAR="
    call :FindNewestShadedJar

    echo The shaded JAR file has been generated at:
    if defined SHADED_JAR (
        echo "!SHADED_JAR!"
    ) else (
        echo metapilot-cli\target\metapilot-cli-*-shaded.jar
    )
    echo.
    echo IMPORTANT: run with "-jar". Do NOT run: java ^<jarfile^>
    echo.
    if defined SHADED_JAR (
        echo To run:
        echo   java -jar "!SHADED_JAR!" --help
    ) else (
        echo To run: java -jar metapilot-cli\target\metapilot-cli-*-shaded.jar --help
    )
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

goto :eof

:FindNewestShadedJar
for /f "delims=" %%F in ('dir /b /a:-d /o:-d "metapilot-cli\target\metapilot-cli-*-shaded.jar" 2^>nul') do (
    set "SHADED_JAR=metapilot-cli\target\%%F"
    goto :eof
)
goto :eof

