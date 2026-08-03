@echo off
setlocal
cd /d "%~dp0"

echo.
echo === NonTunPerc — Validate against recordings ===
echo Working directory: %CD%
echo.

REM Default: open GUI (auto-group from filenames; searches subfolders).
REM Accepts .wav / .aif / .aiff / .flac. Metadata-only grouping (never audio fit).
REM CLI auto-group:
REM   run_validate.bat --cli "C:\path\to\Samples"
REM Manual single-instrument:
REM   run_validate.bat --cli "C:\path\to\folder" cymbal_46cm_medium

if "%~1"=="" (
  echo Mode: GUI
  python validate_against_recordings.py --gui
  exit /b %ERRORLEVEL%
)

if /I "%~1"=="--gui" (
  python validate_against_recordings.py --gui
  exit /b %ERRORLEVEL%
)

if /I "%~1"=="--cli" (
  set "WAV_DIR=%~2"
  set "INSTRUMENT=%~3"
  set "OUT_DIR=%~4"
) else (
  REM Backward compatible: first arg is sample folder
  set "WAV_DIR=%~1"
  set "INSTRUMENT=%~2"
  set "OUT_DIR=%~3"
)

if "%WAV_DIR%"=="" (
  echo [FAIL] Sample folder required for CLI mode.
  echo   run_validate.bat --cli "D:\Samples"
  pause
  exit /b 1
)
if "%OUT_DIR%"=="" set "OUT_DIR=%~dp0validation_out"

if not exist "%WAV_DIR%" (
  echo [FAIL] Sample folder not found: "%WAV_DIR%"
  pause
  exit /b 1
)

echo Mode: CLI
echo Sample dir: %WAV_DIR%
echo Recursive : yes (default)
if "%INSTRUMENT%"=="" (
  echo Grouping  : auto-group metadata-only
) else (
  echo Instrument: %INSTRUMENT% (manual; --no-auto-group)
)
echo Report dir: %OUT_DIR%
echo.

if "%INSTRUMENT%"=="" (
  python validate_against_recordings.py --cli --auto-group --wav-dir "%WAV_DIR%" --out "%OUT_DIR%"
) else (
  python validate_against_recordings.py --cli --no-auto-group --wav-dir "%WAV_DIR%" --instrument %INSTRUMENT% --out "%OUT_DIR%"
)
set ERR=%ERRORLEVEL%

echo.
if %ERR% neq 0 (
  echo [FAIL] validation exited with code %ERR%
) else (
  echo [OK] See "%OUT_DIR%\validation_report.md"
)
echo.
pause
exit /b %ERR%
