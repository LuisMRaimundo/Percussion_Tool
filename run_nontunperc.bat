@echo off
setlocal
cd /d "%~dp0"

echo.
echo === NonTunPerc — non-tuned percussion density model ===
echo Working directory: %CD%
echo.

REM Default: open the GUI. Pass --cli for the headless full pipeline.
if "%~1"=="--cli" (
  echo Mode: headless CLI
  python nontunperc.py --cli %2 %3 %4 %5 %6 %7 %8 %9
) else if "%~1"=="" (
  echo Mode: GUI
  python nontunperc.py
) else (
  python nontunperc.py %*
)

set ERR=%ERRORLEVEL%
echo.
if %ERR% neq 0 (
  echo [FAIL] NonTunPerc exited with code %ERR%
) else (
  echo [OK] NonTunPerc finished.
)
echo.
if not "%~1"=="" pause
exit /b %ERR%
