@echo off
REM Compile and run FindMissingProjects.java
REM Finds projects from projects.txt that don't have .raw files

set JAVA_HOME=C:\Users\daz25hak\AppData\Local\Programs\Eclipse Adoptium\jdk-21.0.8.9-hotspot
set PATH=%JAVA_HOME%\bin;%PATH%

echo.
echo ========================================
echo FindMissingProjects - Find Projects Without .raw Files
echo ========================================
echo.

REM Change to script directory
cd /d "%~dp0"

echo Compiling FindMissingProjects.java...
javac FindMissingProjects.java

if not exist "FindMissingProjects.class" (
    echo.
    echo COMPILATION FAILED!
    echo.
    pause
    exit /b 1
)

echo Compilation successful!
echo.
echo Running FindMissingProjects...
echo ========================================
echo.

java FindMissingProjects

if %ERRORLEVEL% NEQ 0 (
    echo.
    echo ========================================
    echo EXECUTION FAILED with exit code %ERRORLEVEL%
    echo ========================================
    pause
    exit /b %ERRORLEVEL%
)

echo.
echo ========================================
echo Done! Check missing_projects.txt for results.
echo ========================================
pause

