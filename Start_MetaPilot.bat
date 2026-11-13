@echo off
REM Starts MetaPilot.exe from the folder this file is in.
REM Files extracted from a downloaded zip carry Windows' "downloaded from the internet" mark,
REM which makes SmartScreen block the unsigned MetaPilot.exe. The mark is removed from the
REM files in this folder before starting, so later starts are not blocked either.

setlocal
set "APP_DIR=%~dp0"
set "APP_EXE=%APP_DIR%MetaPilot.exe"

if not exist "%APP_EXE%" (
    echo MetaPilot.exe was not found in:
    echo   %APP_DIR%
    echo Keep this file in the same folder as MetaPilot.exe.
    pause
    exit /b 1
)

powershell -NoProfile -ExecutionPolicy Bypass -Command "Get-ChildItem -LiteralPath $env:APP_DIR -Recurse -File -ErrorAction SilentlyContinue | Unblock-File -ErrorAction SilentlyContinue" >nul 2>nul

start "" /D "%APP_DIR%" "%APP_EXE%" %*
endlocal
