@echo off
setlocal enabledelayedexpansion

REM Compile and run ExtractRawFiles.java
REM Extracts .raw file names from PRIDE FTP/HTTP directory links

set JAVA_HOME=C:\Users\daz25hak\AppData\Local\Programs\Eclipse Adoptium\jdk-21.0.8.9-hotspot
set PATH=%JAVA_HOME%\bin;%PATH%

echo.
echo ========================================
echo ExtractRawFiles - PRIDE .raw File Extractor
echo ========================================
echo.

REM Change to script directory
cd /d "%~dp0"

REM Find Apache Commons Net JAR (try Maven local repository)
set COMMONS_NET_JAR=
set FOUND_JAR=0

REM Try specific version first
if exist "%USERPROFILE%\.m2\repository\commons-net\commons-net\3.9.0\commons-net-3.9.0.jar" (
    set COMMONS_NET_JAR=%USERPROFILE%\.m2\repository\commons-net\commons-net\3.9.0\commons-net-3.9.0.jar
    set FOUND_JAR=1
    goto :found_jar
)

REM Try to find any version (search in subdirectories)
for /f "delims=" %%d in ('dir /b /ad "%USERPROFILE%\.m2\repository\commons-net\commons-net" 2^>nul') do (
    if exist "%USERPROFILE%\.m2\repository\commons-net\commons-net\%%d\commons-net-*.jar" (
        for %%f in ("%USERPROFILE%\.m2\repository\commons-net\commons-net\%%d\commons-net-*.jar") do (
            set COMMONS_NET_JAR=%%f
            set FOUND_JAR=1
            goto :found_jar
        )
    )
)

:found_jar
if !FOUND_JAR!==0 (
    echo Warning: Apache Commons Net JAR not found in Maven repository.
    echo FTP functionality will be disabled, but HTTP should still work.
    echo.
    set CLASSPATH=.
) else (
    echo Found Apache Commons Net: !COMMONS_NET_JAR!
    set CLASSPATH=!COMMONS_NET_JAR!
)

echo.
echo Compiling ExtractRawFiles.java...
javac -cp "%CLASSPATH%" ExtractRawFiles.java

REM Check if .class file was created (this is the reliable check)
if not exist "ExtractRawFiles.class" (
    echo.
    echo ========================================
    echo COMPILATION FAILED!
    echo ========================================
    echo.
    echo The ExtractRawFiles.class file was not created.
    echo Please check the compilation errors above.
    echo.
    pause
    exit /b 1
)

echo Compilation successful!
echo.
echo Running ExtractRawFiles...
echo ========================================
echo.

java -cp "%CLASSPATH%;." ExtractRawFiles

set EXIT_CODE=%ERRORLEVEL%

if %EXIT_CODE% NEQ 0 (
    echo.
    echo ========================================
    echo EXECUTION FAILED with exit code %EXIT_CODE%
    echo ========================================
    pause
    exit /b %EXIT_CODE%
)

echo.
echo ========================================
echo Done! Check raw_files_list.txt for results.
echo ========================================
pause
