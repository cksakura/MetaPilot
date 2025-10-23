@echo off
echo Testing MetaPilot logging...
echo.

echo Running MetaPilot for 10 seconds to generate logs...
timeout /t 2 /nobreak > nul

java -jar target\metapilot-1.0-SNAPSHOT.jar --cli --help > nul 2>&1

echo.
echo Checking for log files...
if exist "logs\metapilot.log" (
    echo Found logs\metapilot.log
    echo Content:
    type logs\metapilot.log
) else (
    echo logs\metapilot.log not found
)

if exist "result\logging.log" (
    echo Found result\logging.log
    echo Content:
    type result\logging.log
) else (
    echo result\logging.log not found
)

echo.
echo Test completed.
