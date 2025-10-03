@echo off
echo Searching for MetaPilot reports...
echo.

REM Search for reports directories
echo Looking for 'reports' directories:
for /r %%i in (reports) do (
    if exist "%%i" (
        echo Found: %%i
        dir "%%i\*.html" /b 2>nul
        echo.
    )
)

echo.
echo Looking for HTML report files:
for /r %%i in (*report*.html) do (
    echo Found: %%i
)

echo.
echo Looking for analysis result directories:
for /r %%i in (*analysis*) do (
    if exist "%%i" (
        echo Found: %%i
    )
)

echo.
echo Looking for result directories:
for /r %%i in (*result*) do (
    if exist "%%i" (
        echo Found: %%i
    )
)

pause
